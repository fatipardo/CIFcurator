# CIFcurator server (demo)

A web page that runs the whole CIFcurator workflow from a keyword, without questions, and gives
back a zip file with every file written during the curation.

## Start it

From the repository root, in the CIFcurator conda environment (`environment.yml`: it has MMseqs2 and MUSCLE):

```bash
conda activate CIFcurator
python server/cifcurator_server.py --open
```

and the page opens at <http://127.0.0.1:8000> (without `--open`, open that address yourself). Stop the
server with Ctrl-C. The server only uses Python's standard library: there is nothing else to install.
It uses the CIFcurator code of this checkout (`src/`), even if another version is installed.

The page is served to this computer only. `--port` changes the port, and `--host 0.0.0.0` makes it
reachable from other computers. There is no limit on jobs per person yet, and the server is not
meant to face the internet on its own: a public server should sit behind a web server with HTTPS
(`deploy/` runs it on a cloud server with Docker and Caddy; see `deploy/README.md`).

On your own computer, you probably want the list of recent jobs and no deletion:

```bash
python server/cifcurator_server.py --open --show-recent-jobs --retention-days 0
```

## Several users

There is no sign-up and no cookie. Every job gets a **job number** such as `K7QF-2M9X-HT4P-WN8D`:
16 random characters (about 79 bits, never sequential), without look-alikes (no 0/O, 1/I/L). It is the
only key to the job: whoever has it can see the progress, the step logs and the zip file, and cancel
the job; nobody else can find it. The page shows it after submitting, with a **Copy** button, and the
link `/?job=<number>` can be bookmarked. **Check a job** on the front page opens a job from its number
(any case, with or without dashes).

- **Recent jobs** (the 20 latest jobs of *everybody*, `GET /api/jobs`) are hidden unless the server is
  started with `--show-recent-jobs`: only use it when you are the only user.
- **Retention**: a job (uploaded files, results, zip file) is deleted 7 days after it ends
  (`--retention-days N`, `0` = never). The page shows the date.
- **Limits per job** (default: none): `--max-entries N` refuses a keyword that matches more than N PDB
  entries, right after the search and before anything is downloaded (the count is in the Step 1 log);
  `--max-csm N` allows at most N computed models (blank on the form = N); `--job-timeout-hours H` stops
  a job still running after H hours. The page tells users about the limits that are set.
- **Email** (optional field on the form): one email when the job ends, finished or failed (not when
  cancelled), with the keyword, the job number, the link and the deletion date; no attachment. The
  address is never in `job.json` or the job's page: it is kept in `server_data/notify/<job number>`
  until the email is sent, then deleted. At most 10 emails per address per day (`--max-emails-per-day`).
  A failed email does not fail the job (it is written to the server log).
- Pages are sent with `Referrer-Policy: no-referrer` and `Cache-Control: no-store`, so the job number
  in the address is not passed to other sites or kept in caches.

Emails go through an SMTP server, normally a transactional email service:

```bash
export CIFCURATOR_SMTP_PASSWORD=...        # or --smtp-password-file FILE (outside the code)
python server/cifcurator_server.py --public-url https://cifcurator.example.org \
    --smtp-host email-smtp.us-west-2.amazonaws.com --smtp-user USER --mail-from 'CIFcurator <noreply@example.org>'
```

`--smtp-port` (default 587, 465 with `--smtp-security ssl`), `--smtp-security starttls|ssl|none`. Without
`--smtp-host`, the email is only written to the server log. To read the emails on your own computer, run
[Mailpit](https://mailpit.axllent.org/) and start the server with
`--smtp-host 127.0.0.1 --smtp-port 1025 --smtp-security none`: the emails show up at <http://localhost:8025>.

## What a job does

| Step | What | Choices made for you |
| --- | --- | --- |
| 1 | `CIFcurator_Download.py`: download the structures that match the keyword (`raw_bank/`), clean them (`clean_bank/`) | plain-text keyword; computed models (AlphaFold DB) only if you tick **Include computed structure models** (optionally at most N of them). A computed model has one protein chain only |
| (after 1) | Leave out structures (`excluded_structures.txt`) | only if you list **structures to leave out**: they are not given to Steps 3-5 but stay in `clean_bank/`. Without a list, every structure is curated |
| 3 | MolID: cluster the chains with MMseqs2 (`clustering_storage/`), name the chains (`standard_MolID_bank/`) | your **input file** if you upload one, otherwise the one generated from the clusters (`input_file.txt`, menu option 4). All concatenations are accepted |
| 4 | `CIFcurator_ChainStandardization.py` (`standard_ChainID_bank/`) | your **reference sequences** JSON if you upload one, otherwise the consensus of all the structures (menu option 2 -> 4). Saved as `StandardSequences.json` |
| 5 | `CIFcurator_ResidueStandardization.py` (`standard_ResidueID_bank/`, `OldResID_NewResID_Map.csv`) | occupancy threshold **80** unless you change it; MUSCLE files in `residue_alignments/` (one per chain; `A_b.fasta` = chain Ab, so chains AB and Ab do not clash on macOS) |

The zip file has the project directory: every `*_bank/` folder above, `input_file.txt`,
`clustering_storage/`, `StandardSequences.json`, the track-changes CSV files, and `logs/` (the output of every
step, `parameters.json` with the keyword, threshold and CIFcurator commit, and `error.txt` if a step
failed). If a step fails, a zip file with what was written so far can still be downloaded. The paths in these files
are relative to the project directory (the steps are given relative paths), so a zip file never shows
the job number or where the server keeps its jobs, and it can be shared without giving access to the job.

Structures to leave out are typed in the box (or uploaded as a text file), one or more per line,
separated by spaces or commas, `#` starting a comment. A name can be a file name of `clean_bank/`
(`pdb_00001abc+01.cif`, `.cif` optional), a PDB ID (`1abc`, any case) or an entry (`pdb_00001abc`,
`AF-P00780-F1`); the last two leave out every assembly of the entry. Names that match no downloaded
structure are listed as a warning in `excluded_structures.txt` and next to the step on the page; if
every structure is left out, the job fails.

## Change and rerun

Downloads are slow, so a job that ended (finished, failed or cancelled, as long as Step 1 finished) can
be changed and run again without downloading the structures again:

1. On the job's page, download the files to edit (`input_file.txt`, `StandardSequences.json`,
   `excluded_structures.txt`; links under the steps, also for a failed job) and click **Change and rerun**.
   The form is filled in with the job's settings, and its job number is in **Advanced Options > Start
   from an earlier job** (you can also type a job number there).
2. Change the structures to leave out, upload an edited input file or reference sequences, or change
   the threshold. Files uploaded to the earlier job are used again unless you untick **Use the ... of
   the earlier job** or upload a new one. The keyword and the computed model options must stay as they
   were (otherwise it is not a rerun: run a new job).
3. Before you submit, the form says where the rerun starts: *"This will rerun from Step N; downloads are
   reused."* It starts from the earliest step whose input changed: the structures to leave out -> the
   exclusion step and Steps 3-5; the input file -> Step 3; the reference sequences -> Step 4; the
   threshold -> Step 5. **Step to start from** (Advanced Options) can choose an earlier step, never a
   later one (it would ignore the change). If nothing changed, choose a step to run it again as it is.

The rerun is a **new job with its own job number**; the earlier job and its results are not changed.
The steps not run again show as *reused*. The rerun shares the downloaded structures (and the output of
the steps not run again) of the earlier job by hard links: no copy and no extra disk space, and it keeps
working if the earlier job is deleted first. It is **deleted when the earlier job is** (7 days after the
earlier job ended, also for a rerun of a rerun), so downloads are never kept longer than the retention
time, as the PDB changes every week. If that is less than 24 hours away, the form warns you; the page and
the email show the date. The rerun's page and the API never show the earlier job's number.

Jobs run one at a time, in the order they were submitted. Everything is kept in `server/server_data/`
(`--data-dir` to change it):

- `jobs/<job number>/`: `job.json`, the uploaded files (`inputs/`), the project (`project/`), the zip file.
  Deleted 7 days after the job ends (`--retention-days`).
- `notify/<job number>`: the email address of a job, until its email is sent; `email_counts.json`:
  emails sent today, per address (hashed).
- `pdb_seqres.fasta`: the sequences of every PDB entry (about 400 MB). The first job downloads it, which
  takes a few minutes; the next jobs reuse it, and the first job after every weekly PDB release
  (Wednesday 00:00 UTC) downloads it again. It is not put in the zip files.

## Without the web page

`server/pipeline.py` is what every job runs, and works on its own:

```bash
python server/pipeline.py "subtilisin carlsberg" --projdir ./SUBT --zip SUBT.zip
python server/pipeline.py "subtilisin carlsberg" --projdir ./SUBT3 --csm --max-models 5
python server/pipeline.py "subtilisin carlsberg" --projdir ./SUBT2 \
    --input-file SUBT/input_file.txt --reference-json SUBT/standard_ChainID_bank/StandardSequences.json --threshold 70
python server/pipeline.py "subtilisin carlsberg" --projdir ./SUBT4 --exclude-file leave_out.txt
python server/pipeline.py "subtilisin carlsberg" --projdir ./SUBT5 --reuse-from ./SUBT --start-step residueid --threshold 60
python server/pipeline.py -h
```

## Web API

| Request | Answer |
| --- | --- |
| `POST /api/jobs` (form: `keyword`, `threshold`, optional `csm=on` and `max_models`, optional `exclude` (names), optional `email`, optional files `input_file`, `reference_json`, `exclude_file`; for a rerun, `parent` (job number), optional `start_step` (exclude, molid, chainid, residueid), `keep_input_file=on`, `keep_reference_json=on`) | the new job (JSON, `number` = job number), or `{"error": ...}` with status 400 |
| `GET /api/config` | `show_recent_jobs`, `retention_days` |
| `GET /api/jobs` | the 20 latest jobs, only with `--show-recent-jobs` (404 otherwise) |
| `GET /api/jobs/<id>` | the job: `state` (queued, running, done, failed, cancelled), `steps`, `expires` (deletion time), `rerun` (start step, download time) and `can_rerun`. `<id>` = job number, any case, dashes optional |
| `GET /api/jobs/<id>/log/<step>` | log of a step (download, exclude, molid, chainid, residueid) |
| `GET /api/jobs/<id>/download` | the zip file |
| `GET /api/jobs/<id>/file/<name>` | `input_file.txt`, `StandardSequences.json` (written by the job, else the one uploaded) or `excluded_structures.txt` of the job, to edit for a rerun (`files` in the job lists those it has) |
| `POST /api/rerun-plan` (the form of `POST /api/jobs` with `parent`) | what the rerun would do, nothing is created: `start_step`, `allowed_steps`, `message`, `warning`, `expires`, or `error` |
| `POST /api/jobs/<id>/cancel` | cancels a queued or running job |

```bash
curl -F keyword="subtilisin carlsberg" -F threshold=80 -F input_file=@input_file.txt http://127.0.0.1:8000/api/jobs
```
