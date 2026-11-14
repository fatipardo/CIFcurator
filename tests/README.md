# CIFcurator tests

Run everything from the repository root:

```bash
python -m unittest discover -s tests -t . -v
```

The tests import `CIFcurator` from this checkout's `src/`, so the package does not need
to be installed. They do need its dependencies (biopython, numpy, ...).

## What is tested

| File | What | Needs |
| --- | --- | --- |
| `test_units.py` | Pure functions: residue-name conversion, gap percentages, alignment scoring, seqres search, CIF cleaning (including the OCT2026 path fix), the CIF writers of Steps 3-5, the clustering and input-file generation of Step 3, project directories, the options of `CIFcurator_Download.py` (downloads mocked) | nothing (offline, a few s); one test uses `mmseqs` if it is on `PATH` |
| `test_pipeline.py` | The full workflow, Steps 1 -> 5, on live PDB data, and `CIFcurator_Download.py` | network, `mmseqs` (MMseqs2) and `muscle` (v5) on `PATH` |

Steps 1 and 2 of the pipeline test are the two parts of Step 1 (download, then clean), which users run with
`CIFcurator_Download.py`. They call the library functions of the script, so that only the five entries below are downloaded.
`TestDownloadScript` runs the script itself, on "subtilisin carlsberg" (about 30 entries) and one computed structure model:
it checks that four long-established entries and the model are downloaded, and that every file is cleaned.

**Step 1** searches `pdb_seqres.fasta` (wwPDB beta archive, extended PDB IDs such as `pdb_00001aw1`) for "triosephosphate isomerase", as in the tutorial.
The number of hits grows as the PDB grows, so the test only checks that five
long-established entries (1ag1, 1aw1, 1btm, 1klg, 2y62) are still found. It downloads
them and checks that they are valid mmCIF.

**Steps 2-5** run on those freshly downloaded files, through the same scripts and menu
choices as the tutorials. The checks look at properties of the results, not exact bytes,
so harmless PDB remediation does not break them:

- Step 2: one file per biological assembly (1aw1 -> 4), each holding only its own chains, no atoms lost.
- Steps 2-5: every output file has the categories and key items its data require (atom_type, chem_comp,
  struct_asym, entity_poly_seq, pdbx_poly_seq_scheme, citation.id, ...), and they agree with the atoms
  after chains are renamed and residues renumbered (`_support.mmcif_problems`).
- Step 3: 1klg shows up as the outlier. The entity names in `data/TIM_input_file.txt` still match the PDB files. Chains are renamed only to the assigned IDs.
  The script gets the project directory, with `pdb_seqres.fasta` as `seqres.txt` in it, and clusters the chains with MMseqs2 (NOV2026):
  only the input structures are clustered, all TIM chains fall in one cluster, and option 4 writes an input file with
  every entity of the tutorial file, one chain ID per copy (polymer copies A, Ab, ...), the same chain for both TIM entity names.
  `TestStep3GeneratedInputFile` runs Step 3 again with that generated file (accepting the ligand concatenations) and checks the output.
- Step 4: the 1aw1 chains score 1.0 against the consensus. 1btm (different organism) scores below 0.6.
- Step 5: the CIF numbering (`auth_seq_id` + `pdbx_PDB_ins_code`) matches `OldResID_NewResID_Map.csv`. Numbering increases along each chain. Insertion codes appear for gappy columns. Copies of the same protein get the same residue at the same number.

Each step is a test class and runs the steps before it on demand, so you can run one step alone:

```bash
python -m unittest tests.test_pipeline.TestStep4ChainID -v
```

If a step fails, the later steps are reported as skipped, and the skip message names the
failing step. Each step's full output is in `logs/` inside the work directory (see
`CIFCURATOR_KEEP_TEST_OUTPUT`).

## pdb_seqres.fasta cache

Step 1 needs `pdb_seqres.fasta` (~410 MB, downloaded gzipped from the wwPDB beta archive). The first run downloads it to `tests/.cache/`,
which git ignores. Later runs reuse it until it is 30 days old, then download it again.
The download goes to a temporary file, so an interrupted download is never mistaken for
a complete one.

## Environment variables

| Variable | Effect |
| --- | --- |
| `CIFCURATOR_SKIP_NETWORK=1` | skip the pipeline tests (offline / quick run) |
| `CIFCURATOR_TEST_CACHE=<dir>` | cache directory for `pdb_seqres.fasta` (default `tests/.cache`) |
| `CIFCURATOR_SEQRES_MAX_AGE_DAYS=<n>` | re-download age (default 30; `0` forces a download) |
| `CIFCURATOR_KEEP_TEST_OUTPUT=1` | keep the pipeline work directory and print its path |
| `CIFCURATOR_SCRIPT_TIMEOUT=<s>` | per-script timeout (default 1800). Stops a script that hangs, e.g. waiting for input |

The pipeline work directory is created under `$TMPDIR` with a name that contains "cif".
That is on purpose: it is the regression test for the `.cif` path bug fixed in OCT2026.

Each interactive script runs in its own empty working directory, because Step 5 writes
`<chain>.fa` and `<chain>.fasta` (an underscore before each lowercase letter, e.g. `A_b.fasta` for chain Ab, NOV2026) to the current directory and reuses an existing
`<chain>.fasta`. (NOV2026: Step 4 aligns with Biopython and writes no alignment files.)

