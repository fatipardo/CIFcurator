#!/usr/bin/env python
# coding: utf-8
#
r"""CIFcurator server (demo): curate a structural ensemble from a web page.

Start it from the repository root, in the CIFcurator conda environment (MMseqs2 and MUSCLE on PATH):

  python server/cifcurator_server.py

then open http://127.0.0.1:8000 in a web browser. Type a keyword: the whole workflow (Steps 1-5,
server/pipeline.py) runs without questions, and a zip file with every file it wrote can be
downloaded at the end. Optionally upload the Step 3 input file (chain names of the entities)
and/or the Step 4 Standard Sequences (StandardSequences.json); the residue renumbering threshold
of Step 5 defaults to 80. Computed structure models (AlphaFold DB) can be included in Step 1, and
a list of downloaded structures can be left out of Steps 3-5 (by default all of them are used).

Several people can use the same server: there is no sign-up and no cookie. Every job gets a random
job number (e.g. K7QF-2M9X-HT4P-WN8D), the only way to find the job again: whoever has it can see
the job's progress, logs and zip file, and cancel it. Nobody can list the jobs of others (the
"Recent jobs" list is only shown with --show-recent-jobs, for a server used by one person). A job
is deleted 7 days after it ends (--retention-days). Optionally, one email is sent when a job ends
(--smtp-host ...; without it, the email is only written to the server log).

A job that ended can be changed and run again without downloading the structures again ("Change
and rerun", or the job number in "Start from an earlier job"): the rerun is a new job, with its own
job number, that starts from the earliest step whose input changed (the exclusion list, the Step 3
input file, the Step 4 reference sequences or the Step 5 threshold), or from an earlier step chosen
in Advanced Options. It shares the downloaded structures of the earlier job by hard links, and it is
deleted when the earlier job is, so that the downloads are never kept longer than retention-days.

Jobs run one at a time, in the order they were submitted. Everything is kept in
server/server_data/ (or --data-dir): jobs/<job number>/ for every job, and pdb_seqres.fasta, the PDB
sequence file (about 400 MB) shared by all of them. It is downloaded by the first job, and again by
the first job after every weekly PDB release (Wednesday 00:00 UTC).

Only Python's standard library is needed for the server itself.
"""

import argparse, email.message, email.parser, email.policy, email.utils, hashlib, json, os, queue, re
import secrets, shutil, signal, smtplib, ssl, subprocess, sys, threading, time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

SERVER_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SERVER_DIR)
import pipeline  # noqa: E402

PIPELINE_SCRIPT = os.path.join(SERVER_DIR, 'pipeline.py')

MAX_UPLOAD_BYTES = 20 * 1024 * 1024
MAX_KEYWORD_LENGTH = 200
MAX_CSM_MODELS = 1000
MAX_EXCLUDED = 10000
STRUCTURE_NAME_RE = re.compile(r'^[A-Za-z0-9_.+-]{1,100}$')
# Job numbers: 16 random characters (about 79 bits) without look-alikes (no 0/O, 1/I/L), shown in
# groups of 4 (K7QF-2M9X-HT4P-WN8D), accepted in any case, with or without dashes or spaces.
JOB_ID_ALPHABET = '23456789ABCDEFGHJKMNPQRSTUVWXYZ'
JOB_ID_LENGTH = 16
JOB_ID_RE = re.compile('^[{0}]{{{1}}}$'.format(JOB_ID_ALPHABET, JOB_ID_LENGTH))
DEFAULT_RETENTION_DAYS = 7
CLEANUP_INTERVAL = 3600  # seconds between two looks for expired jobs
MAX_EMAIL_LENGTH = 254
EMAIL_RE = re.compile(r'^[A-Za-z0-9.!#$%&\'*+/=?^_`{|}~-]{1,64}@[A-Za-z0-9](?:[A-Za-z0-9-]*[A-Za-z0-9])?'
                      r'(?:\.[A-Za-z0-9](?:[A-Za-z0-9-]*[A-Za-z0-9])?)+$')
DEFAULT_MAX_EMAILS_PER_DAY = 10
RERUN_WARNING_SECONDS = 86400  # warn when a rerun would be deleted in less than this
# What each step a rerun can start from is called in the messages
STEP_NAMES = {'exclude': 'the step that leaves structures out (then Steps 3-5)', 'molid': 'Step 3',
              'chainid': 'Step 4', 'residueid': 'Step 5'}
CHANGE_NAMES = {'exclude': 'the structures to leave out', 'molid': 'the input file for chain names',
                'chainid': 'the reference sequences', 'residueid': 'the occupancy threshold'}
# Files of a job that users may want to edit and upload again in a rerun (also when the job failed)
# (the first file of the list that exists: the one a step wrote, else the one uploaded)
EDITABLE_FILES = {'input_file.txt': ['project/input_file.txt', 'inputs/input_file.txt'],
                  'StandardSequences.json': ['project/standard_ChainID_bank/StandardSequences.json',
                                             'inputs/StandardSequences.json'],
                  'excluded_structures.txt': ['project/excluded_structures.txt']}


class InputError(ValueError):
    pass


def new_job_id():
    return ''.join(secrets.choice(JOB_ID_ALPHABET) for _ in range(JOB_ID_LENGTH))


def normalize_job_id(text):
    """The job number as stored (K7QF2M9XHT4PWN8D) from what a user typed; KeyError if it cannot be one."""
    job_id = re.sub(r'[\s-]', '', str(text or '')).upper()
    if not JOB_ID_RE.match(job_id):
        raise KeyError(text)
    return job_id


def format_job_id(job_id):
    """The job number as shown: K7QF-2M9X-HT4P-WN8D."""
    return '-'.join(job_id[i:i + 4] for i in range(0, len(job_id), 4))


def fmt_utc(t):
    return time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime(t))


# -- Email ------------------------------------------------------------------------------------
class Mailer:
    """
    Sends the email at the end of a job through an SMTP server (a transactional email service).
    Without an SMTP host, the email is only written to the server log. Never raises: a failed
    email must not fail the job. At most max_per_day emails per address per day (UTC); the
    counts are kept in a file, with the addresses hashed.
    """

    def __init__(self, host=None, port=587, security='starttls', user=None, password=None,
                 sender='CIFcurator <noreply@localhost>', max_per_day=DEFAULT_MAX_EMAILS_PER_DAY,
                 counts_file=None):
        self.host, self.port, self.security = host, port, security
        self.user, self.password, self.sender = user, password, sender
        self.max_per_day = max_per_day
        self.counts_file = counts_file
        self.lock = threading.Lock()

    def count_one(self, address):
        """Counts one email to address today; False if the daily cap is reached."""
        today = time.strftime('%Y-%m-%d', time.gmtime())
        key = hashlib.sha256(address.strip().lower().encode()).hexdigest()
        with self.lock:
            counts = {}
            if self.counts_file:
                try:
                    with open(self.counts_file) as f:
                        counts = json.load(f)
                except (OSError, ValueError):
                    counts = {}
            if counts.get('day') != today:
                counts = {'day': today, 'sent': {}}
            if counts['sent'].get(key, 0) >= self.max_per_day:
                return False
            counts['sent'][key] = counts['sent'].get(key, 0) + 1
            if self.counts_file:
                with open(self.counts_file + '.tmp', 'w') as f:
                    json.dump(counts, f)
                os.replace(self.counts_file + '.tmp', self.counts_file)
            return True

    def send(self, to, subject, body):
        """Returns 'sent', 'logged' (no SMTP host), 'capped' or 'failed'."""
        try:
            if not self.count_one(to):
                log('Email not sent: {0} emails to this address today already.'.format(self.max_per_day))
                return 'capped'
            msg = email.message.EmailMessage()
            msg['From'] = self.sender
            msg['To'] = to
            msg['Subject'] = subject
            msg['Date'] = email.utils.formatdate()
            msg['Message-ID'] = email.utils.make_msgid(domain=email.utils.parseaddr(self.sender)[1].partition('@')[2] or None)
            msg.set_content(body)
            if not self.host:
                log('Email (no --smtp-host, so it is only written here):\n' + msg.as_string())
                return 'logged'
            context = ssl.create_default_context()
            if self.security == 'ssl':
                smtp = smtplib.SMTP_SSL(self.host, self.port, timeout=30, context=context)
            else:
                smtp = smtplib.SMTP(self.host, self.port, timeout=30)
            with smtp:
                if self.security == 'starttls':
                    smtp.starttls(context=context)
                if self.user:
                    smtp.login(self.user, self.password or '')
                smtp.send_message(msg)
            return 'sent'
        except Exception as err:
            log('Email not sent: {0}: {1}'.format(type(err).__name__, err))
            return 'failed'


def log(message):
    print('[{0}] {1}'.format(time.strftime('%Y-%m-%d %H:%M:%S'), message), file=sys.stderr, flush=True)


def job_email(job, public_url, retention_days):
    """Subject and text of the email at the end of a job. Only the keyword comes from the user."""
    keyword = ''.join(c for c in job['keyword'] if c.isprintable())
    if len(keyword) > 60:
        keyword = keyword[:57] + '...'
    number = format_job_id(job['id'])
    done = job['state'] == 'done'
    subject = 'CIFcurator job {0} {1}: {2}'.format(number, 'finished' if done else 'failed', keyword)
    lines = ['Your CIFcurator job has finished.' if done else
             'Your CIFcurator job failed. The page of the job shows which step failed and its log.',
             '',
             'Keyword: ' + keyword,
             'Job number: ' + number,
             '',
             ('Results: ' if done else 'Details: ') + '{0}/?job={1}'.format(public_url.rstrip('/'), number),
             '']
    if job.get('expires') and job.get('rerun'):
        lines.append('The results are kept until {0}, when the earlier job whose downloads this rerun reuses is deleted; '
                     'then they are deleted. Anyone with the job number can view them, so keep it private.'.format(
                         fmt_utc(job['expires'])))
    elif job.get('expires'):
        lines.append('The results are kept for {0} days, until {1}; then they are deleted. Anyone with the job '
                     'number can view them, so keep it private.'.format(retention_days, fmt_utc(job['expires'])))
    else:
        lines.append('Anyone with the job number can view the results, so keep it private.')
    lines += ['',
              'You get this email because this address was given when the job was submitted. It is the only '
              'email about this job, and the address has now been deleted from the server.',
              '',
              '-- ',
              'CIFcurator']
    return subject, '\n'.join(lines) + '\n'


# -- Jobs -------------------------------------------------------------------------------------
class JobStore:
    """
    One directory per job, named by its job number: job.json (parameters and state), inputs/
    (uploaded files), project/ (the CIFcurator project directory), pipeline.log, and the zip file at
    the end. The email address for the end-of-job email is NOT in it: it is kept in notify/<job
    number> until the email is sent, so that it is never shown to whoever has the job number.
    """

    def __init__(self, data_dir, retention_days=DEFAULT_RETENTION_DAYS, mailer=None, public_url='',
                 max_entries=None, job_timeout_hours=None):
        self.data_dir = os.path.abspath(data_dir)
        self.jobs_dir = os.path.join(self.data_dir, 'jobs')
        self.notify_dir = os.path.join(self.data_dir, 'notify')
        self.seqfile = os.path.join(self.data_dir, 'pdb_seqres.fasta')
        os.makedirs(self.jobs_dir, exist_ok=True)
        os.makedirs(self.notify_dir, mode=0o700, exist_ok=True)
        self.retention_days = retention_days  # 0: jobs are never deleted
        self.max_entries = max_entries  # largest number of PDB entries a job may download (None: no limit)
        self.job_timeout = job_timeout_hours * 3600 if job_timeout_hours else None  # seconds (None: no limit)
        self.mailer = mailer or Mailer()
        self.public_url = public_url
        self.lock = threading.Lock()
        self.queue = queue.Queue()
        self.current = None  # (job_id, Popen) of the running job
        self.cancelled = set()
        # Jobs left unfinished by a previous server are not resumed
        for job in self.list():
            if job['state'] in ('queued', 'running'):
                self.finish(job['id'], 'failed', 'The server was stopped before this job finished. Submit it again.')

    def path(self, job_id, *parts):
        return os.path.join(self.jobs_dir, normalize_job_id(job_id), *parts)

    def notify_path(self, job_id):
        return os.path.join(self.notify_dir, normalize_job_id(job_id))

    def get(self, job_id):
        try:
            with open(self.path(job_id, 'job.json')) as f:
                return json.load(f)
        except (OSError, ValueError):
            raise KeyError(job_id)

    def update(self, job_id, **changes):
        with self.lock:
            job = self.get(job_id)
            job.update(changes)
            path = self.path(job_id, 'job.json')
            with open(path + '.tmp', 'w') as f:
                json.dump(job, f, indent=2)
            os.replace(path + '.tmp', path)
            return job

    def list(self):
        jobs = []
        for job_id in os.listdir(self.jobs_dir):
            try:
                jobs.append(self.get(job_id))
            except KeyError:
                pass
        return sorted(jobs, key=lambda j: j['submitted'], reverse=True)

    def read_input(self, job_id, name):
        """Text of an uploaded file of a job (inputs/<name>), or None."""
        try:
            with open(self.path(job_id, 'inputs', name)) as f:
                return f.read()
        except OSError:
            return None

    def download_time(self, job_id):
        """When the structures of a job were downloaded (for a rerun: by the job it reuses)."""
        job = self.get(job_id)
        if job.get('rerun'):
            return job['rerun']['downloaded']
        try:
            with open(self.path(job_id, 'project', 'logs', 'steps.json')) as f:
                return {s['key']: s for s in json.load(f)['steps']}['download'].get('finished')
        except (OSError, ValueError, KeyError):
            return None

    def rerun_plan(self, parent, keyword, threshold, input_file=None, reference_json=None, csm=False, max_models=None,
                   exclude=None, start_step=None, keep_input_file=False, keep_reference_json=False, now=None):
        """
        What a rerun of the job parent (a job number as typed) with these settings would do. Returns a
        dict: start_step, the steps it may start from (allowed_steps: the automatic one, from what
        changed, or an earlier one), the message to show, expires (when it would be deleted: with
        parent) and warning (if that is soon); or 'error' when it cannot be a rerun. input_file and
        reference_json are what was uploaded (None: nothing); with keep_*, None means the file of
        parent. The dict also has the files the rerun uses (input_file, reference_json).
        """
        now = now or time.time()
        try:
            parent_id = normalize_job_id(parent)
            job = self.get(parent_id)
        except KeyError:
            return {'error': 'There is no job with the number {0} to start from.{1}'.format(
                str(parent).strip()[:40], ' Jobs are deleted {0} days after they end.'.format(self.retention_days)
                if self.retention_days else '')}
        if job['state'] in ('queued', 'running'):
            return {'error': 'The earlier job has not ended yet. Wait until it ends to change and rerun it.'}
        if pipeline.step_states(self.path(parent_id, 'project')).get('download') not in ('done', 'reused'):
            return {'error': 'The earlier job did not finish downloading the structures, so there is nothing to reuse. '
                             'Clear "Start from an earlier job" to run a new job.'}
        if job.get('expires') and job['expires'] <= now:
            return {'error': 'The earlier job has expired and is being deleted. Clear "Start from an earlier job" to run a new job.'}
        if keyword != job['keyword'] or bool(csm) != bool(job.get('csm')) or (csm and max_models != job.get('max_models')):
            return {'error': 'A rerun reuses the structures downloaded by the earlier job, so the keyword and the computed '
                             'model options must stay as they were (keyword: {0!r}; computed models: {1}). To change them, '
                             'clear "Start from an earlier job" and run a new job.'.format(
                                 job['keyword'], ('up to {0}'.format(job['max_models']) if job.get('max_models') else 'all')
                                 if job.get('csm') else 'no')}
        old_input, old_reference = self.read_input(parent_id, 'input_file.txt'), self.read_input(parent_id, 'StandardSequences.json')
        if input_file is None and keep_input_file:
            input_file = old_input
        if reference_json is None and keep_reference_json:
            reference_json = old_reference
        changed = []
        if sorted(n.lower() for n in exclude or []) != sorted(n.lower() for n in job.get('exclude') or []):
            changed.append('exclude')
        if input_file != old_input:
            changed.append('molid')
        if reference_json != old_reference:
            changed.append('chainid')
        if threshold != job['threshold']:
            changed.append('residueid')
        auto = changed[0] if changed else None
        allowed = pipeline.RERUN_STEPS[:pipeline.RERUN_STEPS.index(auto) + 1] if auto else list(pipeline.RERUN_STEPS)
        plan = {'parent': parent_id, 'auto_step': auto, 'allowed_steps': allowed, 'changed': changed,
                'expires': job.get('expires'), 'downloaded': self.download_time(parent_id),
                'parent_has_input_file': old_input is not None, 'parent_has_reference_json': old_reference is not None,
                'input_file': input_file, 'reference_json': reference_json}
        if start_step and start_step not in allowed:
            plan['error'] = ('You changed {0}, which {1} uses: the rerun must start from {1} or earlier, otherwise the '
                             'change would be ignored.'.format(CHANGE_NAMES[auto], STEP_NAMES[auto]))
            return plan
        step = start_step or auto
        if step is None:
            plan['error'] = ('Nothing is different from the earlier job yet. Change an input, or choose a step to start '
                             'from in Advanced Options to run it again as it is.')
            return plan
        plan['start_step'] = step
        plan['message'] = 'This will rerun from {0}; downloads are reused.'.format(STEP_NAMES[step])
        if job.get('expires'):
            plan['message'] += (' Like the earlier job, its results will be deleted on {0}.'.format(fmt_utc(job['expires'])))
            if job['expires'] - now < RERUN_WARNING_SECONDS:
                plan['warning'] = ("This job's results will be deleted on {0}. Download them before then.".format(
                    fmt_utc(job['expires'])))
        return plan

    def create(self, keyword, threshold, input_file=None, reference_json=None, csm=False, max_models=None,
               exclude=None, email=None, rerun=None):
        """A new job. rerun: the plan of rerun_plan (without 'error') to start from a step of an earlier job."""
        while True:
            job_id = new_job_id()
            try:
                os.makedirs(self.path(job_id))
                break
            except FileExistsError:  # 1 chance in 10^23, but never reuse a number
                continue
        os.makedirs(self.path(job_id, 'inputs'))
        if rerun:
            # Hard links of the earlier job's downloads (and of the output of the other steps not run again).
            # Under the lock, so that the cleaner cannot delete the earlier job meanwhile.
            with self.lock:
                source = self.path(rerun['parent'], 'project')
                if not os.path.isdir(os.path.join(source, 'clean_bank')):
                    shutil.rmtree(self.path(job_id), ignore_errors=True)
                    raise InputError('The earlier job has just been deleted. Clear "Start from an earlier job" to run a new job.')
                pipeline.reuse_outputs(source, self.path(job_id, 'project'), rerun['start_step'])
        submitted = time.time()
        # The zip file name has no part of the job number: the zip file may be shared
        job = {'id': job_id, 'number': format_job_id(job_id), 'keyword': keyword, 'threshold': threshold,
               'state': 'queued', 'csm': csm, 'max_models': max_models, 'exclude': exclude or [],
               'submitted': submitted, 'input_file': None, 'reference_json': None, 'notify': bool(email),
               'zip_name': 'CIFcurator_{0}_{1}.zip'.format(slug(keyword), time.strftime('%Y%m%d-%H%M%S', time.gmtime(submitted)))}
        if rerun:
            # The earlier job's number is kept (it is never shown: see status) and its expiry is inherited
            job.update(parent=rerun['parent'], expires=rerun['expires'], expires_at_most=rerun['expires'],
                       rerun={'start_step': rerun['start_step'], 'start_title': STEP_NAMES[rerun['start_step']],
                              'downloaded': rerun['downloaded'], 'changed': rerun['changed']})
        if email:
            fd = os.open(self.notify_path(job_id), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            with os.fdopen(fd, 'w') as f:
                f.write(email)
        if input_file is not None:
            with open(self.path(job_id, 'inputs', 'input_file.txt'), 'w') as f:
                f.write(input_file)
            job['input_file'] = 'input_file.txt'
        if reference_json is not None:
            with open(self.path(job_id, 'inputs', 'StandardSequences.json'), 'w') as f:
                f.write(reference_json)
            job['reference_json'] = 'StandardSequences.json'
        if exclude:
            with open(self.path(job_id, 'inputs', 'excluded_structures.txt'), 'w') as f:
                f.write('\n'.join(exclude) + '\n')
        with open(self.path(job_id, 'job.json'), 'w') as f:
            json.dump(job, f, indent=2)
        self.queue.put(job_id)
        return job

    def cancel(self, job_id):
        job_id = normalize_job_id(job_id)
        job = self.get(job_id)
        if job['state'] not in ('queued', 'running'):
            return job
        with self.lock:
            self.cancelled.add(job_id)
            current = self.current
        if current and current[0] == job_id:
            try:
                os.killpg(current[1].pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        return self.finish(job_id, 'cancelled', 'Cancelled by the user.')

    def finish(self, job_id, state, error=None):
        """A job ended (done, failed or cancelled): set its expiry date, send the email if one was asked for."""
        now = time.time()
        job = self.get(job_id)
        expires = now + self.retention_days * 86400 if self.retention_days else None
        if job.get('expires_at_most') is not None:  # a rerun: deleted with the job whose downloads it reuses
            expires = min(expires, job['expires_at_most']) if expires else job['expires_at_most']
        changes = {'state': state, 'finished': now, 'expires': expires}
        if error is not None:
            changes['error'] = error
        job = self.update(job_id, **changes)
        self.notify(job)
        return job

    def notify(self, job):
        """Sends the end-of-job email (done or failed, not cancelled), then deletes the address."""
        path = self.notify_path(job['id'])
        try:
            with open(path) as f:
                address = f.read().strip()
        except OSError:
            return
        try:
            if address and job['state'] in ('done', 'failed'):
                subject, body = job_email(job, self.public_url, self.retention_days)
                self.mailer.send(address, subject, body)
        finally:
            try:
                os.remove(path)
            except OSError:
                pass

    # -- deletion of old jobs ----------------------------------------------------------------
    def delete_expired(self, now=None):
        """Deletes the jobs that ended more than retention_days ago. Returns their job numbers."""
        now = now or time.time()
        deleted = []
        for job in self.list():
            if job['state'] in ('queued', 'running') or not job.get('expires') or job['expires'] > now:
                continue
            with self.lock:  # not while a rerun links its files
                shutil.rmtree(self.path(job['id']), ignore_errors=True)
            deleted.append(job['id'])
        for name in os.listdir(self.notify_dir):  # addresses of jobs that no longer exist
            if not os.path.isdir(os.path.join(self.jobs_dir, name)):
                try:
                    os.remove(os.path.join(self.notify_dir, name))
                except OSError:
                    pass
        if deleted:
            log('Deleted {0} job(s) older than {1} days.'.format(len(deleted), self.retention_days))
        return deleted

    def cleaner(self):
        while True:
            try:
                self.delete_expired()
            except Exception as err:  # never let the cleaner die
                log('Could not delete old jobs: {0}'.format(err))
            time.sleep(CLEANUP_INTERVAL)

    def editable_file(self, job_id, name):
        """Path of a file of EDITABLE_FILES of the job, or None."""
        for rel in EDITABLE_FILES[name]:
            path = self.path(job_id, *rel.split('/'))
            if os.path.isfile(path):
                return path
        return None

    def queue_position(self, job_id):
        waiting = [j for j in self.list() if j['state'] == 'queued']
        waiting.sort(key=lambda j: j['submitted'])
        ids = [j['id'] for j in waiting]
        return ids.index(job_id) + 1 if job_id in ids else None

    # -- worker ------------------------------------------------------------------------------
    def worker(self):
        while True:
            job_id = self.queue.get()
            with self.lock:
                if job_id in self.cancelled:
                    continue
            try:
                self.run(job_id)
            except Exception as err:  # never let the worker die
                self.finish(job_id, 'failed', 'Server error: {0}'.format(err))

    def run(self, job_id):
        job = self.update(job_id, state='running', started=time.time())
        cmd = [sys.executable, PIPELINE_SCRIPT, job['keyword'],
               '--projdir', self.path(job_id, 'project'), '--seqfile', self.seqfile,
               '--threshold', str(job['threshold']), '--zip', self.path(job_id, job['zip_name'])]
        if job.get('csm'):
            cmd.append('--csm')
            if job.get('max_models'):
                cmd += ['--max-models', str(job['max_models'])]
        if job['input_file']:
            cmd += ['--input-file', self.path(job_id, 'inputs', job['input_file'])]
        if job['reference_json']:
            cmd += ['--reference-json', self.path(job_id, 'inputs', job['reference_json'])]
        if job.get('exclude'):
            cmd += ['--exclude-file', self.path(job_id, 'inputs', 'excluded_structures.txt')]
        if job.get('rerun'):
            cmd += ['--start-step', job['rerun']['start_step']]
        if self.max_entries:
            cmd += ['--max-entries', str(self.max_entries)]
        timed_out = False
        with open(self.path(job_id, 'pipeline.log'), 'w') as out:
            proc = subprocess.Popen(cmd, cwd=self.path(job_id), stdout=out, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, start_new_session=True)
            with self.lock:
                self.current = (job_id, proc)
            try:
                proc.wait(timeout=self.job_timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                proc.wait()
            with self.lock:
                self.current = None
        if job_id in self.cancelled:
            return
        if timed_out:
            hours = self.job_timeout / 3600
            self.finish(job_id, 'failed', 'The job was stopped after {0:g} hours, the longest a job may run on this '
                        'server. Please use a more specific keyword, or leave some structures out.'.format(hours))
            return
        if proc.returncode == 0:
            self.finish(job_id, 'done')
        else:
            error = pipeline.tail(self.path(job_id, 'project', 'logs', 'error.txt'), 40) \
                or pipeline.tail(self.path(job_id, 'pipeline.log'), 40)
            self.finish(job_id, 'failed', error.strip())

    def status(self, job_id):
        """What the web page shows: the job and its steps."""
        job_id = normalize_job_id(job_id)
        job = self.get(job_id)
        # The number of the earlier job of a rerun is never shown: this job's number may be shared
        for key in ('parent', 'expires_at_most'):
            job.pop(key, None)
        try:
            with open(self.path(job_id, 'project', 'logs', 'steps.json')) as f:
                steps = json.load(f)['steps']
        except (OSError, ValueError):
            steps = [{'key': k, 'title': t, 'state': 'pending'} for k, t in pipeline.STEPS]
        job['steps'] = steps
        job['queue_position'] = self.queue_position(job_id) if job['state'] == 'queued' else None
        job['zip_ready'] = os.path.isfile(self.path(job_id, job['zip_name']))
        job['seqfile_stale'] = pipeline.seqfile_is_stale(self.seqfile)
        job['files'] = [name for name in EDITABLE_FILES if self.editable_file(job_id, name)]
        job['can_rerun'] = (job['state'] not in ('queued', 'running') and
                            {s['key']: s['state'] for s in steps}.get('download') in ('done', 'reused') and
                            not (job.get('expires') and job['expires'] <= time.time()))
        return job


def slug(text):
    return re.sub(r'[^A-Za-z0-9]+', '_', text).strip('_')[:40] or 'job'


# -- Validation of the form -------------------------------------------------------------------
def clean_keyword(value):
    keyword = ' '.join((value or '').split())
    if not keyword:
        raise InputError('Please type a keyword, e.g. triosephosphate isomerase.')
    if len(keyword) > MAX_KEYWORD_LENGTH:
        raise InputError('The keyword is too long (at most {0} characters).'.format(MAX_KEYWORD_LENGTH))
    return keyword


def clean_threshold(value):
    try:
        threshold = int(str(value).strip() or pipeline.DEFAULT_THRESHOLD)
    except ValueError:
        raise InputError('The occupancy threshold must be a whole number between 1 and 100.')
    if not 1 <= threshold <= 100:
        raise InputError('The occupancy threshold must be between 1 and 100.')
    return threshold


def clean_max_models(value, limit=None):
    """
    Largest number of computed models; blank means all of them (None), or the server's limit
    (--max-csm) if it has one.
    """
    value = str(value or '').strip()
    blank = 'the server limit, {0}'.format(limit) if limit else 'all of them'
    if not value:
        return limit
    try:
        max_models = int(value)
    except ValueError:
        raise InputError('The number of computed models must be a whole number (or blank for {0}).'.format(blank))
    if not 1 <= max_models <= (limit or MAX_CSM_MODELS):
        raise InputError('The number of computed models must be between 1 and {0} (or blank for {1}).'.format(
            limit or MAX_CSM_MODELS, blank))
    return max_models


def clean_exclude(text):
    """Structures to leave out of Steps 3-5: names separated by new lines, commas or spaces."""
    names = pipeline.read_name_list(text)
    bad = [n for n in names if not STRUCTURE_NAME_RE.match(n)]
    if bad:
        raise InputError('These are not names of structure files (letters, digits and _ . + - only, e.g. 1abc, '
                         'pdb_00001abc+01.cif or AF-P00780-F1+01.cif): {0}'.format(', '.join(repr(n[:40]) for n in bad[:5])))
    if len(names) > MAX_EXCLUDED:
        raise InputError('At most {0} structures can be left out.'.format(MAX_EXCLUDED))
    return names


def clean_email(value):
    """Optional address for the end-of-job email; None if blank."""
    address = (value or '').strip()
    if not address:
        return None
    if len(address) > MAX_EMAIL_LENGTH or not EMAIL_RE.match(address):
        raise InputError('This is not an email address: {0!r}. Leave the box empty for no email.'.format(address[:60]))
    return address


def decode_upload(data, what):
    try:
        return data.decode('utf-8-sig')
    except UnicodeDecodeError:
        raise InputError('The {0} is not a text file (UTF-8).'.format(what))


def clean_input_file(data):
    """Step 3 input file: ENTITY NAME:chainID,chainID,... per line. Blank lines and Windows line ends are removed."""
    lines = []
    for number, line in enumerate(decode_upload(data, 'input file').splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        name, sep, chains = line.partition(':')
        if not sep or not name.strip() or not all(c.strip() for c in chains.split(',')):
            raise InputError('Line {0} of the input file is not in the format ENTITY NAME:chainID,chainID,... '
                             '(e.g. TRIOSEPHOSPHATE ISOMERASE:A,Ab): {1!r}'.format(number, line[:100]))
        lines.append(line)
    if not lines:
        raise InputError('The input file is empty.')
    return '\n'.join(lines) + '\n'


def clean_reference_json(data):
    """Step 4 Standard Sequences: a JSON object {chain ID: sequence}."""
    try:
        sequences = json.loads(decode_upload(data, 'reference sequences file'))
    except ValueError as err:
        raise InputError('The reference sequences file is not valid JSON ({0}). It should look like '
                         '{{"A": "MAPSRKFF...", "B": "..."}}, as StandardSequences.json written by Step 4.'.format(err))
    if not isinstance(sequences, dict) or not sequences:
        raise InputError('The reference sequences file must be a JSON object {"chain ID": "sequence", ...}.')
    for chid, seq in sequences.items():
        if not isinstance(seq, str) or not seq.strip() or not chid.strip():
            raise InputError('In the reference sequences file, chain {0!r} needs a sequence (text).'.format(chid))
    return json.dumps(sequences, indent=4) + '\n'


def parse_multipart(content_type, body):
    """Returns {field name: (filename or None, bytes)} of a multipart/form-data body."""
    message = email.parser.BytesParser(policy=email.policy.HTTP).parsebytes(
        b'Content-Type: ' + content_type.encode('latin-1') + b'\r\n\r\n' + body)
    if not message.is_multipart():
        raise InputError('Expected a form (multipart/form-data).')
    fields = {}
    for part in message.iter_parts():
        name = part.get_param('name', header='content-disposition')
        if name:
            fields[name] = (part.get_filename(), part.get_payload(decode=True) or b'')
    return fields


# -- HTTP -------------------------------------------------------------------------------------
class Handler(BaseHTTPRequestHandler):
    store = None  # set in main()
    show_recent_jobs = False  # --show-recent-jobs: GET /api/jobs lists everybody's jobs
    max_csm = None  # --max-csm: largest number of computed models per job (None: MAX_CSM_MODELS)
    server_version = 'CIFcuratorServer/0.2'

    def end_headers(self):
        # The job number is in the page address (?job=...): never send it to other sites, never cache
        self.send_header('Referrer-Policy', 'no-referrer')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('X-Frame-Options', 'DENY')
        self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self' 'unsafe-inline'; "
                         "style-src 'self' 'unsafe-inline'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'")
        self.send_header('Cache-Control', 'no-store')
        super().end_headers()

    def log_message(self, format, *args):
        if not self.path.startswith('/api/jobs/') or self.command != 'GET':
            super().log_message(format, *args)

    def send_json(self, data, status=HTTPStatus.OK):
        body = json.dumps(data).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_error_json(self, status, message):
        self.send_json({'error': message}, status)

    def send_file(self, path, content_type, download_name=None):
        size = os.path.getsize(path)
        self.send_response(HTTPStatus.OK)
        self.send_header('Content-Type', content_type)
        self.send_header('Content-Length', str(size))
        if download_name:
            self.send_header('Content-Disposition', 'attachment; filename="{0}"'.format(download_name))
        self.end_headers()
        with open(path, 'rb') as f:
            shutil.copyfileobj(f, self.wfile)

    def route(self):
        parts = [p for p in self.path.split('?')[0].split('/') if p]
        return parts

    def do_GET(self):
        parts = self.route()
        try:
            if parts in ([], ['index.html']):
                return self.send_file(os.path.join(SERVER_DIR, 'static', 'index.html'), 'text/html; charset=utf-8')
            if parts == ['api', 'config']:
                return self.send_json({'show_recent_jobs': self.show_recent_jobs,
                                       'retention_days': self.store.retention_days,
                                       'max_entries': self.store.max_entries, 'max_csm': self.max_csm,
                                       'job_timeout_hours': self.store.job_timeout and self.store.job_timeout / 3600})
            if parts == ['api', 'jobs'] and self.show_recent_jobs:
                return self.send_json({'jobs': [self.store.status(j['id']) for j in self.store.list()[:20]]})
            if len(parts) == 3 and parts[:2] == ['api', 'jobs']:
                return self.send_json(self.store.status(parts[2]))
            if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] == 'download':
                job = self.store.get(parts[2])
                path = self.store.path(parts[2], job['zip_name'])
                if not os.path.isfile(path):
                    return self.send_error_json(HTTPStatus.NOT_FOUND, 'The zip file is not ready yet.')
                return self.send_file(path, 'application/zip', job['zip_name'])
            if len(parts) == 5 and parts[:2] == ['api', 'jobs'] and parts[3] == 'file':
                if parts[4] not in EDITABLE_FILES:
                    raise KeyError(parts[4])
                path = self.store.editable_file(parts[2], parts[4])
                if not path:
                    return self.send_error_json(HTTPStatus.NOT_FOUND, 'This job has no {0}.'.format(parts[4]))
                return self.send_file(path, 'application/octet-stream', parts[4])
            if len(parts) == 5 and parts[:2] == ['api', 'jobs'] and parts[3] == 'log':
                key = parts[4]
                if key not in dict(pipeline.STEPS):
                    raise KeyError(key)
                path = self.store.path(parts[2], 'project', 'logs', key + '.log')
                if not os.path.isfile(path):
                    return self.send_error_json(HTTPStatus.NOT_FOUND, 'This step has no log yet.')
                return self.send_file(path, 'text/plain; charset=utf-8')
        except KeyError:
            return self.send_error_json(HTTPStatus.NOT_FOUND, 'No such job.')
        self.send_error_json(HTTPStatus.NOT_FOUND, 'Not found.')

    def do_POST(self):
        parts = self.route()
        try:
            if parts == ['api', 'jobs']:
                return self.submit()
            if parts == ['api', 'rerun-plan']:
                return self.submit(dry_run=True)
            if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] == 'cancel':
                return self.send_json(self.store.cancel(parts[2]))
        except KeyError:
            return self.send_error_json(HTTPStatus.NOT_FOUND, 'No such job.')
        except InputError as err:
            return self.send_error_json(HTTPStatus.BAD_REQUEST, str(err))
        self.send_error_json(HTTPStatus.NOT_FOUND, 'Not found.')

    def submit(self, dry_run=False):
        """A new job from the form; with dry_run, only the rerun plan (nothing is created)."""
        length = int(self.headers.get('Content-Length') or 0)
        if length > MAX_UPLOAD_BYTES:
            raise InputError('The uploaded files are too large (at most {0} MB).'.format(MAX_UPLOAD_BYTES // 2**20))
        fields = parse_multipart(self.headers.get('Content-Type', ''), self.rfile.read(length))
        text = lambda name: fields.get(name, (None, b''))[1].decode('utf-8', 'replace')
        keyword = clean_keyword(text('keyword'))
        threshold = clean_threshold(text('threshold'))
        csm = text('csm').strip().lower() in ('on', 'true', '1', 'yes')
        max_models = clean_max_models(text('max_models'), self.max_csm) if csm else None
        input_file = reference_json = None
        filename, data = fields.get('input_file', (None, b''))
        if filename and data:
            input_file = clean_input_file(data)
        filename, data = fields.get('reference_json', (None, b''))
        if filename and data:
            reference_json = clean_reference_json(data)
        exclude_text = text('exclude')
        filename, data = fields.get('exclude_file', (None, b''))
        if filename and data:
            exclude_text += '\n' + decode_upload(data, 'list of structures to leave out')
        exclude = clean_exclude(exclude_text)
        email = clean_email(text('email'))
        parent = text('parent').strip()
        rerun = None
        if parent:
            start_step = text('start_step').strip() or None
            if start_step and start_step not in pipeline.RERUN_STEPS:
                raise InputError('Unknown step to start from: {0!r}.'.format(start_step[:40]))
            rerun = self.store.rerun_plan(parent, keyword, threshold, input_file, reference_json, csm, max_models, exclude,
                                          start_step, text('keep_input_file') == 'on', text('keep_reference_json') == 'on')
            input_file, reference_json = rerun.pop('input_file', input_file), rerun.pop('reference_json', reference_json)
        if dry_run:
            if rerun is None:
                return self.send_json({'error': None, 'message': None})
            plan = dict(rerun)
            plan.pop('parent', None)
            return self.send_json(plan)
        if rerun and rerun.get('error'):
            raise InputError(rerun['error'])
        job = self.store.create(keyword, threshold, input_file, reference_json, csm, max_models, exclude, email, rerun)
        self.send_json(self.store.status(job['id']), HTTPStatus.CREATED)


def check_tools():
    missing = [tool for tool in ('mmseqs', 'muscle') if shutil.which(tool) is None]
    if missing:
        print('WARNING: {0} not found on PATH: Step 3 needs mmseqs (MMseqs2) and Step 5 needs muscle (v5). '
              'Jobs will fail at those steps.'.format(' and '.join(missing)), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                     epilog=__doc__.split('\n\n', 1)[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--host', default='127.0.0.1',
                        help='Address to listen on (default: %(default)s, i.e. only this computer)')
    parser.add_argument('--port', type=int, default=8000, help='Port (default: %(default)s)')
    parser.add_argument('--data-dir', default=os.path.join(SERVER_DIR, 'server_data'),
                        help='Where jobs and the PDB sequence file are kept (default: server/server_data)')
    parser.add_argument('--open', action='store_true', help='Open the page in the web browser')
    parser.add_argument('--show-recent-jobs', action='store_true',
                        help='Show the list of the 20 latest jobs, of EVERYBODY, on the page (GET /api/jobs). '
                             'Only for a server used by one person, e.g. on your own computer. Off by default: '
                             'a job can then only be found with its job number')
    parser.add_argument('--retention-days', type=int, default=DEFAULT_RETENTION_DAYS,
                        help='Jobs (inputs, results, zip file) are deleted this many days after they end '
                             '(default: %(default)s; 0 = never)')
    limits = parser.add_argument_group('limits per job, for a server used by many people (default: none)')
    limits.add_argument('--max-entries', type=int, metavar='N',
                        help='A job whose keyword matches more than N PDB entries stops before downloading anything')
    limits.add_argument('--max-csm', type=int, metavar='N',
                        help='At most N computed models per job (blank on the form = N), instead of {0}'.format(MAX_CSM_MODELS))
    limits.add_argument('--job-timeout-hours', type=float, metavar='H',
                        help='A job still running after H hours is stopped (and fails)')
    parser.add_argument('--public-url',
                        help='Address of the page for the links in the emails, e.g. https://cifcurator.example.org '
                             '(default: the address the server listens on)')
    mail = parser.add_argument_group('email at the end of a job (optional; without --smtp-host the email is '
                                     'only written to the server log)')
    mail.add_argument('--smtp-host', help='SMTP server of the transactional email service')
    mail.add_argument('--smtp-port', type=int, help='SMTP port (default: 587, or 465 with --smtp-security ssl)')
    mail.add_argument('--smtp-security', choices=('starttls', 'ssl', 'none'), default='starttls',
                      help='Encryption (default: %(default)s; none only for a local test server such as Mailpit)')
    mail.add_argument('--smtp-user', help='SMTP user name')
    mail.add_argument('--smtp-password-file',
                      help='File with the SMTP password (default: the environment variable CIFCURATOR_SMTP_PASSWORD). '
                           'Keep it outside the code')
    mail.add_argument('--mail-from', default='CIFcurator <noreply@localhost>',
                      help='Sender of the emails (default: %(default)s)')
    mail.add_argument('--max-emails-per-day', type=int, default=DEFAULT_MAX_EMAILS_PER_DAY,
                      help='Largest number of emails to one address per day (default: %(default)s)')
    args = parser.parse_args()
    if args.retention_days < 0:
        parser.error('--retention-days must be 0 (never delete) or more')
    for option in ('max_entries', 'max_csm', 'job_timeout_hours'):
        if getattr(args, option) is not None and getattr(args, option) <= 0:
            parser.error('--{0} must be more than 0'.format(option.replace('_', '-')))
    if args.max_csm is not None and args.max_csm > MAX_CSM_MODELS:
        parser.error('--max-csm must be at most {0}'.format(MAX_CSM_MODELS))

    password = os.environ.get('CIFCURATOR_SMTP_PASSWORD')
    if args.smtp_password_file:
        with open(args.smtp_password_file) as f:
            password = f.read().strip()
    check_tools()
    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    url = 'http://{0}:{1}'.format('127.0.0.1' if args.host in ('0.0.0.0', '') else args.host, httpd.server_port)
    mailer = Mailer(args.smtp_host, args.smtp_port or (465 if args.smtp_security == 'ssl' else 587),
                    args.smtp_security, args.smtp_user, password, args.mail_from, args.max_emails_per_day,
                    os.path.join(os.path.abspath(args.data_dir), 'email_counts.json'))
    Handler.store = JobStore(args.data_dir, args.retention_days, mailer, args.public_url or url,
                             args.max_entries, args.job_timeout_hours)
    Handler.show_recent_jobs = args.show_recent_jobs
    Handler.max_csm = args.max_csm
    threading.Thread(target=Handler.store.worker, daemon=True).start()
    if args.retention_days:
        threading.Thread(target=Handler.store.cleaner, daemon=True).start()
    print('CIFcurator server running at {0}  (Ctrl-C to stop)'.format(url))
    print('Jobs are kept in {0}, {1}'.format(Handler.store.jobs_dir, 'for {0} days after they end'.format(
        args.retention_days) if args.retention_days else 'until you delete them'))
    print('Limits per job: {0}'.format(', '.join(text for value, text in (
        (args.max_entries, 'at most {0} PDB entries'.format(args.max_entries)),
        (args.max_csm, 'at most {0} computed models'.format(args.max_csm)),
        (args.job_timeout_hours, 'stopped after {0:g} hours'.format(args.job_timeout_hours or 0))) if value) or 'none'))
    print('Recent jobs list: {0}'.format('shown (everybody\'s jobs!)' if args.show_recent_jobs else 'hidden'))
    print('Emails: {0}'.format('sent through {0}:{1}'.format(args.smtp_host, mailer.port) if args.smtp_host
                               else 'only written to this log (no --smtp-host)'))
    if args.open:
        import webbrowser
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print('\nStopping.')
        with Handler.store.lock:
            current = Handler.store.current
        if current:
            try:
                os.killpg(current[1].pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


if __name__ == '__main__':
    main()
