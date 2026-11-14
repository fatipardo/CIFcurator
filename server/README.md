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
reachable from other computers: there is no login and no limit on jobs yet, so only do that on a
trusted network.

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
failed). If a step fails, a zip file with what was written so far can still be downloaded.

Structures to leave out are typed in the box (or uploaded as a text file), one or more per line,
separated by spaces or commas, `#` starting a comment. A name can be a file name of `clean_bank/`
(`pdb_00001abc+01.cif`, `.cif` optional), a PDB ID (`1abc`, any case) or an entry (`pdb_00001abc`,
`AF-P00780-F1`); the last two leave out every assembly of the entry. Names that match no downloaded
structure are listed as a warning in `excluded_structures.txt` and next to the step on the page; if
every structure is left out, the job fails.

To change the result, download the zip file, edit `input_file.txt` (chain names) or
`standard_ChainID_bank/StandardSequences.json` (reference sequences), and upload them with a new job.

Jobs run one at a time, in the order they were submitted. Everything is kept in `server/server_data/`
(`--data-dir` to change it, delete what you no longer need):

- `jobs/<job id>/`: `job.json`, the uploaded files (`inputs/`), the project (`project/`), the zip file.
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
python server/pipeline.py -h
```

## Web API

| Request | Answer |
| --- | --- |
| `POST /api/jobs` (form: `keyword`, `threshold`, optional `csm=on` and `max_models`, optional `exclude` (names), optional files `input_file`, `reference_json`, `exclude_file`) | the new job (JSON), or `{"error": ...}` with status 400 |
| `GET /api/jobs` | the 20 latest jobs |
| `GET /api/jobs/<id>` | the job: `state` (queued, running, done, failed, cancelled), `steps`, end of the running step's log |
| `GET /api/jobs/<id>/log/<step>` | log of a step (download, exclude, molid, chainid, residueid) |
| `GET /api/jobs/<id>/download` | the zip file |
| `POST /api/jobs/<id>/cancel` | cancels a queued or running job |

```bash
curl -F keyword="subtilisin carlsberg" -F threshold=80 -F input_file=@input_file.txt http://127.0.0.1:8000/api/jobs
```
