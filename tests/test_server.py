"""
Tests of the CIFcurator server (server/cifcurator_server.py) used by several people: job numbers,
privacy of the jobs, deletion after the retention time, the email at the end of a job, reruns
of a job from a later step (reusing its downloads), and the limits per job (--max-entries,
--max-csm, --job-timeout-hours).

The pipeline is replaced by a small fake script, and smtplib.SMTP by a fake that records the
messages, so these tests are offline and take a few seconds. TestHTTP starts the server on
127.0.0.1 (a free port) for the duration of the test.
"""

import io
import json
import os
import re
import shutil
import sys
import tempfile
import textwrap
import threading
import time
import unittest
import urllib.error
import urllib.request
from unittest import mock

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'server'))
import cifcurator_server as srv  # noqa: E402

ADDRESS = 'someone@example.org'

FAKE_PIPELINE = textwrap.dedent('''
    import json, os, sys, time, zipfile
    args = sys.argv[1:]
    keyword = args[0]
    if keyword == 'slow':
        time.sleep(60)
    projdir = args[args.index('--projdir') + 1]
    zip_path = args[args.index('--zip') + 1]
    start = args[args.index('--start-step') + 1] if '--start-step' in args else None
    os.makedirs(os.path.join(projdir, 'logs'), exist_ok=True)
    with open(os.path.join(projdir, 'logs', 'args.json'), 'w') as f:
        json.dump(args, f)
    if start is None and keyword != 'nodownload':
        for bank in ('raw_bank', 'clean_bank'):
            os.makedirs(os.path.join(projdir, bank))
            with open(os.path.join(projdir, bank, 'pdb_00001abc+01.cif'), 'w') as f:
                f.write('data_1ABC\\n')
        with open(os.path.join(projdir, 'logs', 'download.log'), 'w') as f:
            f.write('downloaded\\n')
    os.makedirs(os.path.join(projdir, 'standard_MolID_bank'), exist_ok=True)
    with open(os.path.join(projdir, 'input_file.txt'), 'w') as f:
        f.write('PROTEIN:A\\n')
    with open(os.path.join(projdir, 'logs', 'steps.json'), 'w') as f:
        download = 'failed' if keyword == 'nodownload' else 'reused' if start else 'done'
        json.dump({'steps': [{'key': 'download', 'title': 'Step 1', 'state': download, 'finished': 1000.0}]}, f)
    if keyword in ('fail', 'nodownload'):
        with open(os.path.join(projdir, 'logs', 'error.txt'), 'w') as f:
            f.write('Step 3 failed: fake error')
        sys.exit(1)
    with zipfile.ZipFile(zip_path, 'w') as z:
        z.writestr('result.txt', keyword)
''')


class FakeSMTP:
    """Records the messages instead of sending them."""
    sent = []
    fail = False

    def __init__(self, host, port, timeout=None, context=None):
        self.host, self.port = host, port

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def starttls(self, context=None):
        pass

    def login(self, user, password):
        self.user = user

    def send_message(self, msg):
        if FakeSMTP.fail:
            raise OSError('connection refused (fake)')
        FakeSMTP.sent.append(msg)


class ServerTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='cifc-server-test-')
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        script = os.path.join(self.tmp, 'fake_pipeline.py')
        with open(script, 'w') as f:
            f.write(FAKE_PIPELINE)
        patcher = mock.patch.object(srv, 'PIPELINE_SCRIPT', script)
        patcher.start()
        self.addCleanup(patcher.stop)
        smtp = mock.patch.object(srv.smtplib, 'SMTP', FakeSMTP)
        smtp.start()
        self.addCleanup(smtp.stop)
        FakeSMTP.sent = []
        FakeSMTP.fail = False
        self.log = io.StringIO()
        stderr = mock.patch.object(srv.sys, 'stderr', self.log)
        stderr.start()
        self.addCleanup(stderr.stop)
        self.mailer = srv.Mailer('smtp.example.org', 587, 'starttls', 'user', 'secret', 'CIFcurator <noreply@example.org>',
                                 max_per_day=3, counts_file=os.path.join(self.tmp, 'email_counts.json'))
        self.store = srv.JobStore(os.path.join(self.tmp, 'data'), 7, self.mailer, 'https://cifcurator.example.org')

    def submit(self, keyword='subtilisin carlsberg', email=None):
        return self.store.create(keyword, 80, email=email)

    def run_job(self, job):
        self.store.queue.get_nowait()
        self.store.run(job['id'])
        return self.store.get(job['id'])


class TestJobNumbers(unittest.TestCase):
    def test_format(self):
        ids = {srv.new_job_id() for _ in range(2000)}
        self.assertEqual(len(ids), 2000)
        for job_id in ids:
            self.assertRegex(job_id, r'^[2-9A-HJKMNP-Z]{16}$')
            self.assertNotRegex(job_id, '[01OIL]')
        self.assertEqual(srv.format_job_id('K7QF2M9XHT4PWN8D'), 'K7QF-2M9X-HT4P-WN8D')

    def test_typed_forms(self):
        for typed in ('K7QF-2M9X-HT4P-WN8D', 'k7qf2m9xht4pwn8d', ' k7qf 2m9x-ht4p wn8d '):
            self.assertEqual(srv.normalize_job_id(typed), 'K7QF2M9XHT4PWN8D')

    def test_rejected(self):
        for typed in ('', 'K7QF-2M9X-HT4P', 'K7QF2M9XHT4PWN8DX', 'K7QF2M9XHT4PWN80', 'K7QF2M9XHT4PWNOD',
                      '../../etc/passwd', 'K7QF2M9XHT4P/N8D', uuid_like(), None):
            with self.assertRaises(KeyError):
                srv.normalize_job_id(typed)


def uuid_like():
    return 'a' * 32


class TestJobs(ServerTestCase):
    def test_create(self):
        job = self.submit()
        self.assertRegex(job['id'], r'^[2-9A-HJKMNP-Z]{16}$')
        self.assertEqual(job['number'], srv.format_job_id(job['id']))
        self.assertNotIn(job['id'][:4], job['zip_name'])
        self.assertTrue(job['zip_name'].startswith('CIFcurator_subtilisin_carlsberg_'))
        status = self.store.status(job['number'].lower())
        self.assertEqual(status['id'], job['id'])
        self.assertEqual(status['queue_position'], 1)

    def test_unknown_job(self):
        with self.assertRaises(KeyError):
            self.store.status('K7QF-2M9X-HT4P-WN8D')
        with self.assertRaises(KeyError):
            self.store.status('../jobs')

    def test_run_done_and_failed(self):
        job = self.run_job(self.submit())
        self.assertEqual(job['state'], 'done')
        self.assertTrue(self.store.status(job['id'])['zip_ready'])
        self.assertAlmostEqual(job['expires'] - job['finished'], 7 * 86400, delta=1)
        job = self.run_job(self.submit('fail'))
        self.assertEqual(job['state'], 'failed')
        self.assertIn('fake error', job['error'])

    def test_cancel_with_typed_number(self):
        job = self.submit()
        self.store.cancel(job['number'].lower().replace('-', ' '))
        self.assertEqual(self.store.get(job['id'])['state'], 'cancelled')
        self.assertIsNone(self.store.status(job['id'])['queue_position'])

    def test_unfinished_jobs_fail_on_restart(self):
        job = self.submit(email=ADDRESS)
        store = srv.JobStore(os.path.join(self.tmp, 'data'), 7, self.mailer, 'https://cifcurator.example.org')
        job = store.get(job['id'])
        self.assertEqual(job['state'], 'failed')
        self.assertIn('server was stopped', job['error'])
        self.assertEqual(len(FakeSMTP.sent), 1)


class TestRetention(ServerTestCase):
    def test_delete_expired(self):
        old = self.run_job(self.submit('old'))
        new = self.run_job(self.submit('new'))
        queued = self.submit('queued')
        self.assertEqual(self.store.delete_expired(now=old['expires'] - 1), [])
        self.store.update(new['id'], expires=old['expires'] + 3600)
        self.assertEqual(self.store.delete_expired(now=old['expires'] + 1), [old['id']])
        self.assertFalse(os.path.exists(self.store.path(old['id'])))
        with self.assertRaises(KeyError):
            self.store.status(old['id'])
        self.assertEqual(self.store.get(new['id'])['state'], 'done')
        # queued and running jobs are never deleted
        self.assertEqual(self.store.delete_expired(now=time.time() + 100 * 86400), [new['id']])
        self.assertEqual(self.store.get(queued['id'])['state'], 'queued')

    def test_never_delete(self):
        store = srv.JobStore(os.path.join(self.tmp, 'data0'), 0, self.mailer)
        job = store.create('x', 80)
        store.queue.get_nowait()
        store.run(job['id'])
        self.assertIsNone(store.get(job['id'])['expires'])
        self.assertEqual(store.delete_expired(now=time.time() + 1000 * 86400), [])

    def test_leftover_address_deleted(self):
        job = self.submit(email=ADDRESS)
        shutil.rmtree(self.store.path(job['id']))
        self.store.delete_expired()
        self.assertEqual(os.listdir(self.store.notify_dir), [])


class TestEmail(ServerTestCase):
    def job_files_text(self, job_id):
        text = ''
        for root, _, files in os.walk(self.store.path(job_id)):
            for name in files:
                with open(os.path.join(root, name), 'rb') as f:
                    text += f.read().decode('utf-8', 'replace')
        return text

    def test_done(self):
        job = self.submit(email=ADDRESS)
        self.assertTrue(job['notify'])
        self.assertNotIn(ADDRESS, json.dumps(self.store.status(job['id'])))
        self.assertNotIn(ADDRESS, self.job_files_text(job['id']))
        self.assertTrue(os.path.isfile(self.store.notify_path(job['id'])))
        job = self.run_job(job)
        self.assertEqual(len(FakeSMTP.sent), 1)
        msg = FakeSMTP.sent[0]
        body = msg.get_content()
        self.assertEqual(msg['To'], ADDRESS)
        self.assertIn('finished', msg['Subject'])
        self.assertIn(job['number'], msg['Subject'])
        self.assertIn('subtilisin carlsberg', msg['Subject'])
        self.assertIn('Job number: ' + job['number'], body)
        self.assertIn('https://cifcurator.example.org/?job=' + job['number'], body)
        self.assertIn('kept for 7 days', body)
        self.assertIn(srv.fmt_utc(job['expires']), body)
        self.assertFalse(os.path.exists(self.store.notify_path(job['id'])))

    def test_failed(self):
        self.run_job(self.submit('fail', email=ADDRESS))
        self.assertEqual(len(FakeSMTP.sent), 1)
        self.assertIn('failed', FakeSMTP.sent[0]['Subject'])
        self.assertNotIn('fake error', FakeSMTP.sent[0].get_content())

    def test_cancelled_no_email(self):
        job = self.submit(email=ADDRESS)
        self.store.cancel(job['id'])
        self.assertEqual(FakeSMTP.sent, [])
        self.assertFalse(os.path.exists(self.store.notify_path(job['id'])))

    def test_no_address(self):
        job = self.run_job(self.submit())
        self.assertFalse(job['notify'])
        self.assertEqual(FakeSMTP.sent, [])

    def test_daily_cap(self):
        for _ in range(5):
            self.run_job(self.submit(email=ADDRESS))
        self.assertEqual(len(FakeSMTP.sent), 3)
        self.assertIn('3 emails to this address today', self.log.getvalue())
        self.run_job(self.submit(email='other@example.org'))
        self.assertEqual(len(FakeSMTP.sent), 4)
        with open(self.mailer.counts_file) as f:
            self.assertNotIn(ADDRESS, f.read())

    def test_send_failure(self):
        FakeSMTP.fail = True
        job = self.run_job(self.submit(email=ADDRESS))
        self.assertEqual(job['state'], 'done')
        self.assertIn('Email not sent', self.log.getvalue())
        self.assertFalse(os.path.exists(self.store.notify_path(job['id'])))

    def test_no_smtp_host(self):
        self.store.mailer = srv.Mailer()
        self.run_job(self.submit(email=ADDRESS))
        self.assertEqual(FakeSMTP.sent, [])
        self.assertIn('only written here', self.log.getvalue())
        self.assertIn('Subject: CIFcurator job', self.log.getvalue())

    def test_keyword_in_subject_is_safe(self):
        subject, body = srv.job_email({'id': 'K7QF2M9XHT4PWN8D', 'keyword': 'x\r\nBcc: a@b.c' + 'y' * 100,
                                       'state': 'done', 'expires': None}, 'http://h', 7)
        self.assertNotIn('\n', subject)
        self.assertTrue(subject.endswith('...'))
        self.assertLessEqual(len(subject.partition(': ')[2]), 60)

    def test_clean_email(self):
        self.assertIsNone(srv.clean_email('  '))
        self.assertEqual(srv.clean_email(' a.b+c@cs.example.edu '), 'a.b+c@cs.example.edu')
        for bad in ('a', 'a@b', 'a b@c.org', 'a@b.org\nBcc: x@y.org', 'a@b.org, c@d.org', '<a@b.org>', 'x' * 250 + '@b.org'):
            with self.assertRaises(srv.InputError):
                srv.clean_email(bad)


class TestHTTP(ServerTestCase):
    def start(self, show_recent_jobs=False):
        handler = type('H', (srv.Handler,), {'store': self.store, 'show_recent_jobs': show_recent_jobs,
                                            'log_message': lambda *a: None})
        self.httpd = srv.ThreadingHTTPServer(('127.0.0.1', 0), handler)
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        self.addCleanup(self.httpd.server_close)
        self.addCleanup(self.httpd.shutdown)
        self.url = 'http://127.0.0.1:{0}'.format(self.httpd.server_port)

    def get(self, path):
        try:
            with urllib.request.urlopen(self.url + path) as res:
                return res.status, res.headers, res.read()
        except urllib.error.HTTPError as err:
            return err.code, err.headers, err.read()

    def post_form(self, fields, path='/api/jobs'):
        boundary = 'testboundary'
        body = ''.join('--{0}\r\nContent-Disposition: form-data; name="{1}"\r\n\r\n{2}\r\n'.format(boundary, k, v)
                       for k, v in fields.items()) + '--{0}--\r\n'.format(boundary)
        req = urllib.request.Request(self.url + path, data=body.encode(), method='POST',
                                     headers={'Content-Type': 'multipart/form-data; boundary=' + boundary})
        try:
            with urllib.request.urlopen(req) as res:
                return res.status, json.loads(res.read())
        except urllib.error.HTTPError as err:
            return err.code, json.loads(err.read())

    def test_no_job_list(self):
        self.start()
        self.submit()
        status, headers, body = self.get('/api/jobs')
        self.assertEqual(status, 404)
        self.assertNotIn(b'subtilisin', body)
        status, _, body = self.get('/api/config')
        self.assertEqual(json.loads(body), {'show_recent_jobs': False, 'retention_days': 7, 'max_entries': None,
                                            'max_csm': None, 'job_timeout_hours': None})

    def test_job_list_when_asked(self):
        self.start(show_recent_jobs=True)
        self.submit()
        status, _, body = self.get('/api/jobs')
        self.assertEqual(status, 200)
        self.assertEqual(len(json.loads(body)['jobs']), 1)
        self.assertTrue(json.loads(self.get('/api/config')[2])['show_recent_jobs'])

    def test_submit_and_check(self):
        self.start()
        status, job = self.post_form({'keyword': 'subtilisin carlsberg', 'threshold': '80', 'email': ADDRESS})
        self.assertEqual(status, 201)
        self.assertTrue(job['notify'])
        self.assertNotIn(ADDRESS, json.dumps(job))
        for typed in (job['number'], job['number'].lower(), job['id']):
            status, headers, body = self.get('/api/jobs/' + typed)
            self.assertEqual(status, 200)
            self.assertEqual(json.loads(body)['id'], job['id'])
            self.assertNotIn(ADDRESS.encode(), body)
            self.assertEqual(headers['Referrer-Policy'], 'no-referrer')
            self.assertEqual(headers['Cache-Control'], 'no-store')
        self.assertEqual(self.get('/api/jobs/K7QF-2M9X-HT4P-WN8D')[0], 404)
        self.assertEqual(self.get('/api/jobs/K7QF-2M9X-HT4P-WN8D/download')[0], 404)
        self.assertEqual(self.get('/api/jobs/..%2F..%2Fetc/download')[0], 404)

    def test_bad_email(self):
        self.start()
        status, data = self.post_form({'keyword': 'x', 'email': 'not an address'})
        self.assertEqual(status, 400)
        self.assertIn('email address', data['error'])

    def test_page(self):
        self.start()
        status, headers, body = self.get('/?job=K7QF-2M9X-HT4P-WN8D')
        self.assertEqual(status, 200)
        page = body.decode()
        for text in ('id="check-form"', 'name="email"', 'id="copy"', "It's the only way to find it again",
                     'id="recent-card" hidden', '/api/config'):
            self.assertIn(text, page)
        self.assertIn("frame-ancestors 'none'", headers['Content-Security-Policy'])


class TestRerun(ServerTestCase):
    INPUT = 'PROTEIN:A,Ab\n'
    REFERENCE = srv.clean_reference_json(b'{"A": "MAPSRKFF"}')

    def parent(self, **kw):
        job = self.store.create(kw.pop('keyword', 'subtilisin carlsberg'), kw.pop('threshold', 80), **kw)
        return self.run_job(job)

    def plan(self, parent, **kw):
        args = dict(keyword='subtilisin carlsberg', threshold=80)
        args.update(kw)
        return self.store.rerun_plan(parent['number'], **args)

    def rerun(self, parent, **kw):
        plan = self.plan(parent, **kw)
        self.assertNotIn('error', plan)
        keep = {k: kw.get(k) for k in ('keyword', 'threshold', 'csm', 'max_models', 'exclude')}
        job = self.store.create(keep['keyword'] or 'subtilisin carlsberg', keep['threshold'] or 80, plan['input_file'],
                                plan['reference_json'], bool(keep['csm']), keep['max_models'], keep['exclude'], rerun=plan)
        return job, plan

    def test_auto_step(self):
        parent = self.parent(exclude=['1xyz'])
        self.assertEqual(self.plan(parent, threshold=70, exclude=['1xyz'])['start_step'], 'residueid')
        self.assertEqual(self.plan(parent, reference_json=self.REFERENCE, exclude=['1XYZ'])['start_step'], 'chainid')
        self.assertEqual(self.plan(parent, input_file=self.INPUT, exclude=['1xyz'])['start_step'], 'molid')
        self.assertEqual(self.plan(parent, exclude=[])['start_step'], 'exclude')
        plan = self.plan(parent, threshold=50, input_file=self.INPUT, exclude=['1xyz'])
        self.assertEqual((plan['start_step'], plan['changed']), ('molid', ['molid', 'residueid']))
        self.assertEqual(plan['allowed_steps'], ['exclude', 'molid'])
        self.assertIn('rerun from Step 3; downloads are reused', plan['message'])

    def test_manual_step(self):
        parent = self.parent()
        plan = self.plan(parent, threshold=70, start_step='exclude')
        self.assertEqual(plan['start_step'], 'exclude')
        plan = self.plan(parent, reference_json=self.REFERENCE, start_step='residueid')
        self.assertIn('must start from Step 4 or earlier', plan['error'])
        self.assertEqual(plan['allowed_steps'], ['exclude', 'molid', 'chainid'])
        # Nothing changed: any step, but one must be chosen
        self.assertIn('Nothing is different', self.plan(parent)['error'])
        self.assertEqual(self.plan(parent, start_step='residueid')['start_step'], 'residueid')

    def test_kept_files(self):
        parent = self.parent(input_file=self.INPUT, reference_json=self.REFERENCE)
        plan = self.plan(parent, threshold=60, keep_input_file=True, keep_reference_json=True)
        self.assertEqual((plan['start_step'], plan['input_file'], plan['reference_json']), ('residueid', self.INPUT, self.REFERENCE))
        self.assertTrue(plan['parent_has_input_file'])
        # Not kept: back to automatic chain names, a change for Step 3
        self.assertEqual(self.plan(parent, threshold=60, keep_reference_json=True)['start_step'], 'molid')

    def test_not_a_rerun(self):
        parent = self.parent()
        self.assertIn('no job with the number', self.store.rerun_plan('K7QF-2M9X-HT4P-WN8D', 'subtilisin carlsberg', 80)['error'])
        self.assertIn('keyword and the computed', self.plan(parent, keyword='other', threshold=70)['error'])
        self.assertIn('keyword and the computed', self.plan(parent, csm=True, threshold=70)['error'])
        queued = self.store.create('subtilisin carlsberg', 80)
        self.assertIn('not ended yet', self.plan(queued, threshold=70)['error'])
        failed = self.parent(keyword='nodownload')
        self.assertEqual(failed['state'], 'failed')
        self.assertIn('did not finish downloading', self.plan(failed, keyword='nodownload', threshold=70)['error'])
        self.store.update(parent['id'], expires=time.time() - 1)
        self.assertIn('expired', self.plan(parent, threshold=70)['error'])

    def test_failed_parent_can_be_rerun(self):
        parent = self.parent(keyword='fail')
        self.assertEqual(parent['state'], 'failed')
        self.assertTrue(self.store.status(parent['id'])['can_rerun'])
        self.assertEqual(self.plan(parent, keyword='fail', input_file=self.INPUT)['start_step'], 'molid')

    def test_rerun_shares_downloads(self):
        parent = self.parent()
        job, plan = self.rerun(parent, threshold=70)
        self.assertNotEqual(job['id'], parent['id'])
        old = self.store.path(parent['id'], 'project', 'clean_bank', 'pdb_00001abc+01.cif')
        new = self.store.path(job['id'], 'project', 'clean_bank', 'pdb_00001abc+01.cif')
        self.assertEqual(os.stat(old).st_ino, os.stat(new).st_ino)
        self.assertTrue(os.path.isfile(self.store.path(job['id'], 'project', 'standard_MolID_bank', '..', 'input_file.txt')))
        self.assertFalse(os.path.exists(self.store.path(job['id'], 'project', 'standard_ResidueID_bank')))
        job = self.run_job(job)
        self.assertEqual(job['state'], 'done')
        with open(self.store.path(job['id'], 'project', 'logs', 'args.json')) as f:
            args = json.load(f)
        self.assertEqual(args[args.index('--start-step') + 1], 'residueid')
        self.assertEqual(args[args.index('--threshold') + 1], '70')
        # The earlier job is unchanged, and its deletion does not touch the rerun
        self.assertEqual(self.store.get(parent['id'])['state'], 'done')
        shutil.rmtree(self.store.path(parent['id']))
        self.assertTrue(os.path.isfile(new))

    def test_parent_number_never_shown(self):
        parent = self.parent()
        job, _ = self.rerun(parent, threshold=70)
        self.assertEqual(self.store.get(job['id'])['parent'], parent['id'])
        status = json.dumps(self.store.status(job['id']))
        self.assertNotIn(parent['id'], status)
        self.assertNotIn(parent['number'], status)
        self.assertEqual(self.store.status(job['id'])['rerun']['start_step'], 'residueid')

    def test_expiry_inherited(self):
        parent = self.parent()
        self.store.update(parent['id'], expires=time.time() + 3 * 86400)
        job, plan = self.rerun(parent, threshold=70)
        self.assertNotIn('warning', plan)
        self.assertEqual(job['expires'], self.store.get(parent['id'])['expires'])
        job = self.run_job(job)
        self.assertAlmostEqual(job['expires'], time.time() + 3 * 86400, delta=5)
        # A rerun of the rerun is deleted with the first job too
        job2, _ = self.rerun(job, threshold=60)
        self.assertEqual(self.run_job(job2)['expires'], job['expires'])

    def test_warning_when_close(self):
        parent = self.parent()
        expires = time.time() + 5 * 3600
        self.store.update(parent['id'], expires=expires)
        plan = self.plan(parent, threshold=70)
        self.assertIn("This job's results will be deleted on " + srv.fmt_utc(expires), plan['warning'])
        self.assertIn('Download them before then', plan['warning'])

    def test_no_retention(self):
        store = srv.JobStore(os.path.join(self.tmp, 'data0'), 0, self.mailer)
        parent = store.create('x', 80)
        store.queue.get_nowait()
        store.run(parent['id'])
        plan = store.rerun_plan(parent['id'], 'x', 70)
        self.assertIsNone(plan['expires'])
        self.assertNotIn('deleted', plan['message'])

    def test_rerun_email(self):
        parent = self.parent()
        plan = self.plan(parent, threshold=70)
        job = self.store.create('subtilisin carlsberg', 70, email=ADDRESS, rerun=plan)
        self.run_job(job)
        body = FakeSMTP.sent[-1].get_content()
        self.assertIn('when the earlier job whose downloads this rerun reuses is deleted', body)
        self.assertNotIn(parent['number'], body)


class TestRerunHTTP(ServerTestCase):
    start, get, post_form = TestHTTP.start, TestHTTP.get, TestHTTP.post_form

    def test_rerun(self):
        self.start()
        parent = self.run_job(self.submit())
        form = {'keyword': 'subtilisin carlsberg', 'threshold': '65', 'parent': parent['number'].lower()}
        status, plan = self.post_form(form, '/api/rerun-plan')
        self.assertEqual((status, plan['start_step']), (200, 'residueid'))
        self.assertNotIn(parent['id'], json.dumps(plan))
        self.assertEqual(self.store.list()[0]['id'], parent['id'])  # nothing created
        status, plan = self.post_form(dict(form, start_step='nonsense'), '/api/rerun-plan')
        self.assertEqual(status, 400)
        status, job = self.post_form(dict(form, keyword='other'))
        self.assertEqual(status, 400)
        self.assertIn('keyword', job['error'])
        status, job = self.post_form(form)
        self.assertEqual(status, 201)
        self.assertEqual(job['rerun']['start_step'], 'residueid')
        self.assertNotIn(parent['id'], json.dumps(job))
        status, plan = self.post_form({'keyword': 'x'}, '/api/rerun-plan')
        self.assertEqual((status, plan['message']), (200, None))

    def test_editable_files(self):
        self.start()
        job = self.run_job(self.submit('fail'))
        status = self.store.status(job['id'])
        self.assertIn('input_file.txt', status['files'])
        code, headers, body = self.get('/api/jobs/{0}/file/input_file.txt'.format(job['number']))
        self.assertEqual((code, body), (200, b'PROTEIN:A\n'))
        self.assertIn('attachment', headers['Content-Disposition'])
        self.assertEqual(self.get('/api/jobs/{0}/file/job.json'.format(job['number']))[0], 404)
        self.assertEqual(self.get('/api/jobs/{0}/file/StandardSequences.json'.format(job['number']))[0], 404)
        # An uploaded file, when the step wrote none
        job = self.store.create('fail', 80, reference_json='{"A": "MK"}\n')
        job = self.run_job(job)
        code, _, body = self.get('/api/jobs/{0}/file/StandardSequences.json'.format(job['number']))
        self.assertEqual((code, body), (200, b'{"A": "MK"}\n'))

    def test_page(self):
        self.start()
        page = self.get('/')[2].decode()
        for text in ('name="parent"', 'name="start_step"', 'id="rerun"', '/api/rerun-plan', 'id="rerun-warning"',
                     'data-for="input_file"', 'data-for="reference_json"', 'data-for="exclude_file"',
                     '[hidden] { display:none !important; }'):
            self.assertIn(text, page)


class TestPipelineRerun(unittest.TestCase):
    """pipeline.py --start-step / --reuse-from, with the steps replaced by fakes."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='cifc-rerun-test-')
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.old = os.path.join(self.tmp, 'old')
        for out in ('raw_bank/a.cif', 'clean_bank/a.cif', 'clean_bank/b.cif', 'excluded_structures.txt',
                    'logs/excluded_files.txt', 'input_file.txt', 'clustering_storage/x/c.txt', 'standard_MolID_bank/a.cif',
                    'standard_ChainID_bank/a.cif', 'standard_ChainID_bank/StandardSequences.json',
                    'standard_ResidueID_bank/a.cif', 'residue_alignments/A.fasta', 'logs/download.log', 'logs/molid.log'):
            path = os.path.join(self.old, out)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, 'w') as f:
                f.write('b.cif\n' if out.endswith('excluded_files.txt') else out)

    def run_from(self, step, exclude=()):
        new = tempfile.mkdtemp(prefix='new-' + step, dir=self.tmp)
        srv.pipeline.reuse_outputs(self.old, new, step)
        p = srv.pipeline.Pipeline('kw', new, os.path.join(self.tmp, 'seq.fasta'), start_step=step, exclude=exclude)
        ran = []
        for key in ('download', 'exclude_structures', 'molid', 'chainid', 'residueid'):
            setattr(p, key, lambda key=key: ran.append(key))
        p.run_all()
        return new, p, ran

    def test_from_step5(self):
        new, p, ran = self.run_from('residueid')
        self.assertEqual(ran, ['residueid'])
        self.assertEqual([p.status[k]['state'] for k in ('download', 'exclude', 'molid', 'chainid')], ['reused'] * 4)
        self.assertEqual(p.excluded_files, ['b.cif'])
        self.assertEqual(os.stat(os.path.join(new, 'standard_ChainID_bank', 'a.cif')).st_ino,
                         os.stat(os.path.join(self.old, 'standard_ChainID_bank', 'a.cif')).st_ino)
        self.assertTrue(os.path.isfile(os.path.join(new, 'clustering_storage', 'x', 'c.txt')))
        for out in ('standard_ResidueID_bank', 'residue_alignments'):
            self.assertFalse(os.path.exists(os.path.join(new, out)))
        with open(os.path.join(new, 'logs', 'parameters.json')) as f:
            self.assertEqual(json.load(f)['start_step'], 'residueid')

    def test_from_step3(self):
        new, p, ran = self.run_from('molid')
        self.assertEqual(ran, ['molid', 'chainid', 'residueid'])
        for out in ('input_file.txt', 'clustering_storage', 'standard_MolID_bank', 'logs/molid.log'):
            self.assertFalse(os.path.exists(os.path.join(new, out)), out)
        self.assertTrue(os.path.isfile(os.path.join(new, 'excluded_structures.txt')))

    def test_from_exclusion(self):
        new, p, ran = self.run_from('exclude', exclude=['1abc'])
        self.assertEqual(ran, ['exclude_structures', 'molid', 'chainid', 'residueid'])
        self.assertFalse(os.path.exists(os.path.join(new, 'excluded_structures.txt')))
        new, p, ran = self.run_from('exclude')  # the list was emptied
        self.assertEqual(p.status['exclude']['state'], 'skipped')
        self.assertEqual(ran, ['molid', 'chainid', 'residueid'])

    def test_missing_downloads(self):
        shutil.rmtree(os.path.join(self.old, 'clean_bank'))
        with self.assertRaises(srv.pipeline.StepFailed):
            self.run_from('molid')


class TestRelativePaths(unittest.TestCase):
    """The steps get paths relative to their working directory: the job number (in the path of projdir)
    must not end up in the files of the zip (logs, alignment headers, records)."""

    NUMBER = 'K7QF2M9XHT4PWN8D'

    def test_no_job_number_in_commands(self):
        tmp = tempfile.mkdtemp(prefix='cifc-relpath-test-')
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        job = os.path.join(tmp, 'jobs', self.NUMBER)
        projdir = os.path.join(job, 'project')
        for path in ('project/clean_bank/a.cif', 'project/standard_ResidueID_bank/a.cif', 'inputs/input_file.txt',
                     'inputs/StandardSequences.json'):
            os.makedirs(os.path.dirname(os.path.join(job, path)), exist_ok=True)
            open(os.path.join(job, path), 'w').close()
        seqfile = os.path.join(tmp, 'pdb_seqres.fasta')
        open(seqfile, 'w').close()
        p = srv.pipeline.Pipeline('kw', projdir, seqfile, input_file=os.path.join(job, 'inputs', 'input_file.txt'),
                                  reference_json=os.path.join(job, 'inputs', 'StandardSequences.json'), exclude=['1abc'])
        os.makedirs(p.logdir)
        p.excluded_files = ['pdb_00001abc+01.cif']
        calls = []
        p.run = lambda key, cmd, cwd, answers=None, append=False: calls.append((cmd, cwd, answers or ''))
        for step in (p.download, p.molid, p.chainid, p.residueid):
            step()
        self.assertEqual(len(calls), 4)
        for cmd, cwd, answers in calls:
            args = cmd[2:]  # not python and the script of this checkout (not in the job's folder)
            self.assertNotIn(self.NUMBER, ' '.join(args) + answers)
            self.assertFalse([a for a in args if os.path.isabs(a)], cmd)
        # Each relative path leads to the right place from the working directory of its step
        download, molid, chainid, residueid = calls
        here = lambda call, arg: os.path.realpath(os.path.join(call[1], arg))
        self.assertEqual(here(download, '.'), os.path.realpath(projdir))
        self.assertEqual(here(download, download[0][download[0].index('--seqfile') + 1]), os.path.realpath(seqfile))
        self.assertEqual(here(molid, molid[0][molid[0].index('--input-file') + 1]),
                         os.path.realpath(os.path.join(job, 'inputs', 'input_file.txt')))
        self.assertTrue(os.path.isfile(here(molid, molid[0][molid[0].index('--exclude-file') + 1])))
        self.assertEqual(here(residueid, residueid[0][-1]), os.path.realpath(os.path.join(projdir, 'standard_ResidueID_bank')))
        self.assertIn('../inputs/StandardSequences.json', calls[2][2])
        self.assertIn('../standard_ChainID_bank', calls[3][0])

        # The command line in the log, and parameters.json
        p2 = srv.pipeline.Pipeline('kw', projdir, seqfile)
        p2.run('chainid', [sys.executable, '-c', 'pass'], cwd=projdir)
        with open(os.path.join(p2.logdir, 'chainid.log')) as f:
            self.assertEqual(f.readline(), '$ python -c pass\n')
        p.run = p2.run
        for key in ('download', 'exclude_structures', 'molid', 'chainid', 'residueid'):
            setattr(p, key, lambda: None)
        p.run_all()
        with open(os.path.join(p.logdir, 'parameters.json')) as f:
            text = f.read()
        self.assertNotIn(self.NUMBER, text)
        self.assertEqual(json.loads(text)['input_file'], '../inputs/input_file.txt')


class TestLimits(ServerTestCase):
    """Limits per job for a server used by many people: PDB entries, computed models, run time."""

    SEQRES = ('>pdb_00001abc_A mol:protein length:8  TRIOSEPHOSPHATE ISOMERASE\nMKVLAAGG\n'
              '>pdb_00002abc_A mol:protein length:8  UNNAMED PROTEIN\nMKVLAAGG\n'  # same sequence: downloaded too
              '>pdb_00003abc_B mol:protein length:8  TRIOSEPHOSPHATE ISOMERASE\nWWWWWWWW\n'
              '>pdb_00004abc_A mol:protein length:8  LYSOZYME\nKVFGRCEL\n')

    def count(self, keyword, max_entries):
        seqfile = os.path.join(self.tmp, 'pdb_seqres.fasta')
        with open(seqfile, 'w') as f:
            f.write(self.SEQRES)
        import subprocess
        return subprocess.run([sys.executable, srv.pipeline.__file__, '_count', keyword, seqfile, '--max-entries',
                               str(max_entries)], env=srv.pipeline.package_env(), capture_output=True, text=True)

    def test_count_entries(self):
        res = self.count('triosephosphate isomerase', 3)
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertIn('matches 3 PDB entries', res.stdout)
        res = self.count('triosephosphate isomerase', 2)
        self.assertEqual(res.returncode, 1)
        self.assertIn('matches 3 PDB entries, more than the 2 this server can curate', res.stderr)
        self.assertEqual(self.count('nothing like this', 1).returncode, 0)

    def test_pipeline_checks_before_download(self):
        seqfile = os.path.join(self.tmp, 'seq.fasta')
        with open(seqfile, 'w') as f:
            f.write(self.SEQRES)
        p = srv.pipeline.Pipeline('kw', os.path.join(self.tmp, 'proj'), seqfile, max_entries=50, update_seqfile=True)
        calls = []
        p.run = lambda key, cmd, cwd, answers=None, append=False: calls.append((cmd, append))
        with self.assertRaises(srv.pipeline.StepFailed):  # nothing downloaded by the fake run
            p.download()
        (count, append1), (download, append2) = calls
        self.assertEqual(count[2:4], ['_count', 'kw'])
        self.assertEqual(count[-3:], ['--max-entries', '50', '--update'])
        self.assertNotIn('--update', download)  # the sequence file was already downloaded by _count
        self.assertEqual((append1, append2), (False, True))
        calls.clear()
        p.max_entries = None
        with self.assertRaises(srv.pipeline.StepFailed):
            p.download()
        self.assertEqual(len(calls), 1)
        self.assertIn('--update', calls[0][0])

    def test_server_passes_max_entries(self):
        job = self.run_job(self.submit())
        with open(self.store.path(job['id'], 'project', 'logs', 'args.json')) as f:
            self.assertNotIn('--max-entries', json.load(f))
        self.store.max_entries = 2000
        job = self.run_job(self.submit())
        with open(self.store.path(job['id'], 'project', 'logs', 'args.json')) as f:
            args = json.load(f)
        self.assertEqual(args[args.index('--max-entries') + 1], '2000')

    def test_timeout(self):
        self.store.job_timeout = 0.5
        start = time.time()
        job = self.run_job(self.submit('slow'))
        self.assertLess(time.time() - start, 20)
        self.assertEqual(job['state'], 'failed')
        self.assertIn('stopped after 0.000138889 hours', job['error'])
        self.assertIsNone(self.store.current)
        self.store.job_timeout = None
        self.assertEqual(self.run_job(self.submit())['state'], 'done')

    def test_out_of_memory_message(self):
        p = srv.pipeline.Pipeline('kw', os.path.join(self.tmp, 'proj'), os.path.join(self.tmp, 'seq.fasta'))
        os.makedirs(p.logdir)
        with self.assertRaisesRegex(srv.pipeline.StepFailed, 'ran out of memory'):
            p.run('chainid', [sys.executable, '-c', 'import os, signal; os.kill(os.getpid(), signal.SIGKILL)'], cwd=p.projdir)

    def test_max_models(self):
        self.assertIsNone(srv.clean_max_models(''))
        self.assertEqual(srv.clean_max_models('', 50), 50)
        self.assertEqual(srv.clean_max_models('20', 50), 20)
        with self.assertRaisesRegex(srv.InputError, 'between 1 and 50 .or blank for the server limit, 50'):
            srv.clean_max_models('51', 50)
        self.assertEqual(srv.clean_max_models('1000'), 1000)


class TestLimitsHTTP(TestHTTP):
    def test_limits(self):
        self.store.max_entries, self.store.job_timeout = 2000, 4 * 3600
        handler = type('H', (srv.Handler,), {'store': self.store, 'max_csm': 50, 'log_message': lambda *a: None})
        self.httpd = srv.ThreadingHTTPServer(('127.0.0.1', 0), handler)
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        self.addCleanup(self.httpd.server_close)
        self.addCleanup(self.httpd.shutdown)
        self.url = 'http://127.0.0.1:{0}'.format(self.httpd.server_port)
        config = json.loads(self.get('/api/config')[2])
        self.assertEqual((config['max_entries'], config['max_csm'], config['job_timeout_hours']), (2000, 50, 4))
        status, job = self.post_form({'keyword': 'subtilisin', 'csm': 'on', 'max_models': ''})
        self.assertEqual((status, job['max_models']), (201, 50))
        status, err = self.post_form({'keyword': 'subtilisin', 'csm': 'on', 'max_models': '200'})
        self.assertEqual(status, 400)
        self.assertIn('between 1 and 50', err['error'])
        self.assertIn('limits-note', self.get('/')[2].decode())


if __name__ == '__main__':
    unittest.main()
