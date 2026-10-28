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
| `test_units.py` | Pure functions: residue-name conversion, gap percentages, alignment scoring, seqres search, CIF cleaning (including the OCT2026 path fix), the CIF writers of Steps 3-5, project directories | nothing (offline, ~1 s) |
| `test_pipeline.py` | The full workflow, Steps 1 -> 5, on live PDB data | network, `muscle` (v5) on `PATH` |

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
| `CIFCURATOR_SCRIPT_TIMEOUT=<s>` | per-script timeout (default 1800). If MUSCLE fails, Steps 4-5 would otherwise wait forever |

The pipeline work directory is created under `$TMPDIR` with a name that contains "cif".
That is on purpose: it is the regression test for the `.cif` path bug fixed in OCT2026.

Each interactive script runs in its own empty working directory, because the alignment
code writes `Seq.fa`, `<chain>.fa` and `<chain>.fasta` to the current directory and reuses
an existing `<chain>.fasta`.

## On Sherlock

Run the tests in a job, not on a login node, for example:

```bash
srun -p dev -t 00:30:00 -c 2 --mem=4G python -m unittest discover -s tests -t . -v
```
