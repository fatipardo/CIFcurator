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

Jobs run one at a time, in the order they were submitted. Everything is kept in
server/server_data/ (or --data-dir): jobs/<job id>/ for every job, and pdb_seqres.fasta, the PDB
sequence file (about 400 MB) shared by all of them. It is downloaded by the first job, and again by
the first job after every weekly PDB release (Wednesday 00:00 UTC).

Only Python's standard library is needed for the server itself.
"""

import argparse, email.parser, email.policy, json, os, queue, re, shutil, signal, subprocess, sys
import threading, time, uuid
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

SERVER_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SERVER_DIR)
import pipeline  # noqa: E402

MAX_UPLOAD_BYTES = 20 * 1024 * 1024
MAX_KEYWORD_LENGTH = 200
MAX_CSM_MODELS = 1000
MAX_EXCLUDED = 10000
STRUCTURE_NAME_RE = re.compile(r'^[A-Za-z0-9_.+-]{1,100}$')
JOB_ID_RE = re.compile(r'^[0-9a-f]{32}$')


class InputError(ValueError):
    pass


# -- Jobs -------------------------------------------------------------------------------------
class JobStore:
    """
    One directory per job: job.json (parameters and state), inputs/ (uploaded files),
    project/ (the CIFcurator project directory), pipeline.log, and the zip file at the end.
    """

    def __init__(self, data_dir):
        self.data_dir = os.path.abspath(data_dir)
        self.jobs_dir = os.path.join(self.data_dir, 'jobs')
        self.seqfile = os.path.join(self.data_dir, 'pdb_seqres.fasta')
        os.makedirs(self.jobs_dir, exist_ok=True)
        self.lock = threading.Lock()
        self.queue = queue.Queue()
        self.current = None  # (job_id, Popen) of the running job
        self.cancelled = set()
        # Jobs left unfinished by a previous server are not resumed
        for job in self.list():
            if job['state'] in ('queued', 'running'):
                self.update(job['id'], state='failed', finished=time.time(),
                            error='The server was stopped before this job finished. Submit it again.')

    def path(self, job_id, *parts):
        if not JOB_ID_RE.match(job_id):
            raise KeyError(job_id)
        return os.path.join(self.jobs_dir, job_id, *parts)

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

    def create(self, keyword, threshold, input_file=None, reference_json=None, csm=False, max_models=None,
               exclude=None):
        job_id = uuid.uuid4().hex
        os.makedirs(self.path(job_id, 'inputs'))
        job = {'id': job_id, 'keyword': keyword, 'threshold': threshold, 'state': 'queued',
               'csm': csm, 'max_models': max_models, 'exclude': exclude or [],
               'submitted': time.time(), 'input_file': None, 'reference_json': None,
               'zip_name': 'CIFcurator_{0}_{1}.zip'.format(slug(keyword), job_id[:8])}
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
        return self.update(job_id, state='cancelled', finished=time.time(), error='Cancelled by the user.')

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
                self.update(job_id, state='failed', finished=time.time(), error='Server error: {0}'.format(err))

    def run(self, job_id):
        job = self.update(job_id, state='running', started=time.time())
        cmd = [sys.executable, os.path.join(SERVER_DIR, 'pipeline.py'), job['keyword'],
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
        with open(self.path(job_id, 'pipeline.log'), 'w') as log:
            proc = subprocess.Popen(cmd, cwd=self.path(job_id), stdout=log, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, start_new_session=True)
            with self.lock:
                self.current = (job_id, proc)
            proc.wait()
            with self.lock:
                self.current = None
        if job_id in self.cancelled:
            return
        if proc.returncode == 0:
            self.update(job_id, state='done', finished=time.time())
        else:
            error = pipeline.tail(self.path(job_id, 'project', 'logs', 'error.txt'), 40) \
                or pipeline.tail(self.path(job_id, 'pipeline.log'), 40)
            self.update(job_id, state='failed', finished=time.time(), error=error.strip())

    def status(self, job_id):
        """What the web page shows: the job and its steps."""
        job = self.get(job_id)
        try:
            with open(self.path(job_id, 'project', 'logs', 'steps.json')) as f:
                steps = json.load(f)['steps']
        except (OSError, ValueError):
            steps = [{'key': k, 'title': t, 'state': 'pending'} for k, t in pipeline.STEPS]
        job['steps'] = steps
        job['queue_position'] = self.queue_position(job_id) if job['state'] == 'queued' else None
        job['zip_ready'] = os.path.isfile(self.path(job_id, job['zip_name']))
        job['seqfile_stale'] = pipeline.seqfile_is_stale(self.seqfile)
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


def clean_max_models(value):
    """Largest number of computed models; blank means all of them (None)."""
    value = str(value or '').strip()
    if not value:
        return None
    try:
        max_models = int(value)
    except ValueError:
        raise InputError('The number of computed models must be a whole number (or blank for all of them).')
    if not 1 <= max_models <= MAX_CSM_MODELS:
        raise InputError('The number of computed models must be between 1 and {0} (or blank for all of them).'.format(MAX_CSM_MODELS))
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
    server_version = 'CIFcuratorServer/0.1'

    def log_message(self, format, *args):
        if not self.path.startswith('/api/jobs/') or self.command != 'GET':
            super().log_message(format, *args)

    def send_json(self, data, status=HTTPStatus.OK):
        body = json.dumps(data).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
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
            if parts == ['api', 'jobs']:
                return self.send_json({'jobs': [self.store.status(j['id']) for j in self.store.list()[:20]]})
            if len(parts) == 3 and parts[:2] == ['api', 'jobs']:
                return self.send_json(self.store.status(parts[2]))
            if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] == 'download':
                job = self.store.get(parts[2])
                path = self.store.path(parts[2], job['zip_name'])
                if not os.path.isfile(path):
                    return self.send_error_json(HTTPStatus.NOT_FOUND, 'The zip file is not ready yet.')
                return self.send_file(path, 'application/zip', job['zip_name'])
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
            if len(parts) == 4 and parts[:2] == ['api', 'jobs'] and parts[3] == 'cancel':
                return self.send_json(self.store.cancel(parts[2]))
        except KeyError:
            return self.send_error_json(HTTPStatus.NOT_FOUND, 'No such job.')
        except InputError as err:
            return self.send_error_json(HTTPStatus.BAD_REQUEST, str(err))
        self.send_error_json(HTTPStatus.NOT_FOUND, 'Not found.')

    def submit(self):
        length = int(self.headers.get('Content-Length') or 0)
        if length > MAX_UPLOAD_BYTES:
            raise InputError('The uploaded files are too large (at most {0} MB).'.format(MAX_UPLOAD_BYTES // 2**20))
        fields = parse_multipart(self.headers.get('Content-Type', ''), self.rfile.read(length))
        text = lambda name: fields.get(name, (None, b''))[1].decode('utf-8', 'replace')
        keyword = clean_keyword(text('keyword'))
        threshold = clean_threshold(text('threshold'))
        csm = text('csm').strip().lower() in ('on', 'true', '1', 'yes')
        max_models = clean_max_models(text('max_models')) if csm else None
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
        job = self.store.create(keyword, threshold, input_file, reference_json, csm, max_models, exclude)
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
    args = parser.parse_args()

    check_tools()
    Handler.store = JobStore(args.data_dir)
    threading.Thread(target=Handler.store.worker, daemon=True).start()
    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    url = 'http://{0}:{1}'.format('127.0.0.1' if args.host in ('0.0.0.0', '') else args.host, httpd.server_port)
    print('CIFcurator server running at {0}  (Ctrl-C to stop)'.format(url))
    print('Jobs are kept in {0}'.format(Handler.store.jobs_dir))
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        