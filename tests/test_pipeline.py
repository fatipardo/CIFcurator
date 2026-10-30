"""
End-to-end test of the CIFcurator workflow, Steps 1 -> 5, on live PDB data.

Step 1 searches the current pdb_seqres.fasta of the wwPDB beta archive as in the Step 1 tutorial. The number of
hits grows over time, so the test does not check the total. It checks that a fixed
set of long-established TIM entries is still found, downloads only those, and runs
Steps 2-5 on the freshly downloaded files. This checks that the code still works
with what the PDB serves today.

Requirements: network access, and `muscle` (v5) on PATH (Steps 4 and 5).
pdb_seqres.fasta (~410 MB) is cached in tests/.cache and re-downloaded after 30 days;
see _support.py for the environment variables that control this.

Run from the repository root:
    python -m unittest tests.test_pipeline -v

Each step is a test class. A step runs the earlier steps it needs (once per run), so
any single class can be run alone, e.g. tests.test_pipeline.TestStep4ChainID. If an
earlier step fails, the later steps are skipped and the skip message names it.
"""
import contextlib
import glob
import json
import os
import shutil
import string
import sys
import tempfile
import unittest
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _support  # noqa: E402  (makes CIFcurator importable from src/)

from CIFcurator import cifcurator_io, cleanutils, pdbutils  # noqa: E402
from CIFcurator import cifcuratormolidcifutils as molidutils  # noqa: E402

KEYWORD = 'triosephosphate isomerase'
# Tutorial subset: four TIM structures plus 1klg, an HLA complex that only carries a
# TIM peptide and must surface as the outlier in Step 3.
# OCT2026: extended PDB IDs, as in the wwPDB beta archive (1aw1 -> pdb_00001aw1)
PDB_IDS = ['pdb_00001ag1', 'pdb_00001aw1', 'pdb_00001btm', 'pdb_00001klg', 'pdb_00002y62']
OUTLIER = 'pdb_00001klg'
# Biological assemblies, one file each from Step 1 on (1aw1 has four).
ASSEMBLIES = ['pdb_00001ag1+01', 'pdb_00001aw1+01', 'pdb_00001aw1+02', 'pdb_00001aw1+03', 'pdb_00001aw1+04',
              'pdb_00001btm+01', 'pdb_00001klg+01', 'pdb_00002y62+01']
# Step 3 assignment of entity name -> new chain IDs, as in the Step 3 tutorials. The 2y62 assembly
# is a dimer generated from a one-chain asymmetric unit: its second TIM chain is Ab (copies of a
# polymer chain A are named Ab, Ac, ...), and its extra ligand copies get M-P.
MOLID_INPUT = os.path.join(_support.DATA_DIR, 'TIM_input_file.txt')
OCCUPANCY_THRESHOLD = 75
# Same alphabet, in the same order, as get_gap_letter() in Step 5.
GAP_ALPHABET = string.ascii_uppercase + string.ascii_lowercase + string.digits


def entry_id(meta):
    # The extended PDB ID of a fasta description line ('>pdb_00001aw1_A mol:protein ...' -> 'pdb_00001aw1')
    return pdbutils.get_idset_from_metadata([meta])[0]


class Pipeline:
    """Runs the steps lazily, in order, each at most once per test run."""
    STEPS = ['step1', 'step2', 'step3', 'step4', 'step5']

    def __init__(self):
        self.root = None
        self.done = set()
        self.failed = {}

    def require(self, step):
        if _support.SKIP_NETWORK:
            raise unittest.SkipTest('CIFCURATOR_SKIP_NETWORK is set')
        for s in self.STEPS[:self.STEPS.index(step) + 1]:
            if s in self.done:
                continue
            if s not in self.failed:
                try:
                    getattr(self, 'run_' + s)()
                    self.done.add(s)
                    continue
                except Exception as e:
                    self.failed[s] = e
                    raise
            if s == step:
                raise self.failed[s]
            raise unittest.SkipTest('{0} failed: {1}'.format(s, str(self.failed[s]).splitlines()[0]))

    def path(self, *parts):
        return os.path.join(self.root, *parts)

    def log(self, name):
        """Context manager that sends the library's printing to a log file."""
        f = open(self.path('logs', name), 'w')
        stack = contextlib.ExitStack()
        stack.enter_context(f)
        stack.enter_context(contextlib.redirect_stdout(f))
        return stack

    def run_script(self, step, script, args, answers):
        """Run an interactive script in its own empty working directory (see _support.run_script)."""
        _support.run_script(script, args, answers, cwd=self.path('cwd_' + step),
                            log_path=self.path('logs', step + '.log'))

    def run_step1(self):
        # "cif" in the directory name on purpose: regression test for the OCT2026 path fix.
        self.root = tempfile.mkdtemp(prefix='cifcurator_specific_')
        os.mkdir(self.path('logs'))
        self.seqfile = _support.seqres_file()
        with self.log('step1.log'):
            ref_sequences, ref_metadata = pdbutils.retrieve_sequence_from_PDB(
                KEYWORD, mode='metadata', update=False, seqfile=self.seqfile)
            self.keyword_ids = set(pdbutils.get_idset_from_metadata(ref_metadata))
            # OCT2026: as in the tutorial, a sequence search for every reference hit, now in a single
            # reading of pdb_seqres.fasta (before: one reading per hit, close to an hour)
            _, meta = pdbutils.retrieve_sequences_matching_from_PDB(ref_sequences, self.seqfile)
            self.sequence_metadata = list(meta)
            self.sequence_ids = set(pdbutils.get_idset_from_metadata(self.sequence_metadata))
            cifcurator_io.check_project(projdir=self.root, level='raw_bank')
            pdbutils.download_pdb_from_metadata(
                [m for m in self.sequence_metadata if entry_id(m) in PDB_IDS], projdir=self.root)

    def run_step2(self):
        with self.log('step2.log'):
            cifcurator_io.check_project(projdir=self.root, level='clean_bank')
            # OCT2026: clean and simplify are a single step
            cleanutils.process(projdir=self.root, step='clean', source='raw_bank', target='clean_bank')

    def run_step3(self):
        # As in the tutorial: drop the outlier, then assign chain IDs from the input file.
        os.mkdir(self.path('clean_bank_sub'))
        for f in glob.glob(self.path('clean_bank', '*.cif')):
            if not os.path.basename(f).startswith(OUTLIER):
                shutil.copy(f, self.path('clean_bank_sub'))
        os.mkdir(self.path('standard_MolID_bank'))
        # 3 = enter input file, A = write track-changes CSV, 7 = continue, 6 = finalize
        self.run_script('step3', 'CIFcurator_MolID.py',
                        [self.path('clean_bank_sub'), self.path('standard_MolID_bank')],
                        '3\n{0}\nA\n7\n6\n'.format(MOLID_INPUT))

    def run_step4(self):
        os.mkdir(self.path('standard_ChainID_bank'))
        # 2 = consensus menu, 4 = generate consensus, 4 = alignment menu, 4 = align and rename
        self.run_script('step4', 'CIFcurator_ChainStandardization.py',
                        [self.path('standard_MolID_bank'), self.path('standard_ChainID_bank')],
                        '2\n4\n4\n4\n')

    def run_step5(self):
        os.mkdir(self.path('standard_ResidueID_bank'))
        # 1 = MSA menu, 4 = align with threshold, 4 = save template, 3 = write CIFs
        self.run_script('step5', 'CIFcurator_ResidueStandardization.py',
                        [self.path('standard_ChainID_bank'), self.path('standard_ResidueID_bank')],
                        '1\n4\n{0}\n4\n3\nQUIT\n'.format(OCCUPANCY_THRESHOLD))


PIPELINE = Pipeline()


def setUpModule():
    _support.quiet_known_warnings()
    if not _support.SKIP_NETWORK and shutil.which('muscle') is None:
        print('\nWARNING: muscle not found on PATH; Steps 4 and 5 will fail.', file=sys.stderr)


def tearDownModule():
    if PIPELINE.root:
        if _support.KEEP_OUTPUT:
            print('\nPipeline output kept in {0}'.format(PIPELINE.root), file=sys.stderr)
        else:
            shutil.rmtree(PIPELINE.root, ignore_errors=True)


def cif_files(bank):
    return sorted(os.path.basename(f)[:-4] for f in glob.glob(PIPELINE.path(bank, '*.cif')))


def residues(path):
    """Ordered list of (chain, label_seq_id, auth_seq_id, ins_code, comp_id) for polymer residues of a CIF."""
    group, asym, seq, auth, ins, comp = _support.cif_columns(
        path, '_atom_site.group_PDB', '_atom_site.label_asym_id', '_atom_site.label_seq_id',
        '_atom_site.auth_seq_id', '_atom_site.pdbx_PDB_ins_code', '_atom_site.label_comp_id')
    out = []
    for row in zip(group, asym, seq, auth, ins, comp):
        if row[0] == 'ATOM' and (not out or out[-1] != row[1:]):
            out.append(row[1:])
    return out


def assert_consistent(test, bank):
    """OCT2026 reviewer: every file must have the categories and key items its data require."""
    for name in cif_files(bank):
        with test.subTest(bank=bank, structure=name):
            test.assertEqual(_support.mmcif_problems(PIPELINE.path(bank, name + '.cif')), [])


def residue_order(num, ins):
    """Sort key for a Step 5 residue number: 45 < 45 A < 45 B < ... < 45 Z < 45 a < ... < 46."""
    letters = '' if ins in ('?', '.', '') else ins
    return (int(num), len(letters), [GAP_ALPHABET.index(c) for c in letters])


class TestStep1Download(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        PIPELINE.require('step1')

    def test_seqres_cache_is_fasta(self):
        with open(PIPELINE.seqfile) as f:
            # OCT2026: the beta archive file uses extended PDB IDs
            self.assertTrue(f.readline().startswith('>pdb_'))
        self.assertGreater(os.path.getsize(PIPELINE.seqfile), 100 * 1024 * 1024)

    def test_keyword_search_finds_known_entries(self):
        missing = set(PDB_IDS) - PIPELINE.keyword_ids
        self.assertFalse(missing, 'not found by keyword {0!r}: {1}'.format(KEYWORD, sorted(missing)))

    def test_sequence_search_finds_known_entries(self):
        missing = set(PDB_IDS) - PIPELINE.sequence_ids
        self.assertFalse(missing, 'not found by sequence search: {0}'.format(sorted(missing)))
        self.assertTrue(all(m.startswith('>') for m in PIPELINE.sequence_metadata))

    def test_downloaded_files_are_valid_mmcif(self):
        # One file per biological assembly, straight from the wwPDB assemblies archive, saved as PDBID+NN.
        names = cif_files('raw_bank')
        self.assertEqual(names, ASSEMBLIES)
        for name in names:
            with self.subTest(name=name):
                (atoms,) = _support.cif_columns(PIPELINE.path('raw_bank', name + '.cif'), '_atom_site.id')
                self.assertGreater(len(atoms), 100)


class TestStep2Clean(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        PIPELINE.require('step2')

    def test_one_file_per_biological_assembly(self):
        # Also the regression test for the OCT2026 path fix: the work directory contains "cif".
        self.assertEqual(cif_files('clean_bank'), ASSEMBLIES)

    def test_assembly_files_keep_every_atom(self):
        # Each wwPDB assembly file already is one assembly: Step 2 rewrites it whole, under the same name.
        for name in ASSEMBLIES:
            with self.subTest(assembly=name):
                raw = _support.cif_columns(PIPELINE.path('raw_bank', name + '.cif'),
                                           '_atom_site.id', '_atom_site.label_asym_id')
                clean = _support.cif_columns(PIPELINE.path('clean_bank', name + '.cif'),
                                             '_atom_site.id', '_atom_site.label_asym_id')
                self.assertGreater(len(clean[0]), 0)
                self.assertEqual(clean, raw)

    def test_pdb_id_restored(self):
        for name in ASSEMBLIES:
            with self.subTest(assembly=name):
                with open(PIPELINE.path('clean_bank', name + '.cif')) as f:
                    text = f.read()
                self.assertIn("_entry.id   '{0}'".format(name.split('+')[0]), text)

    def test_consistent_with_mmcif_dictionary(self):
        assert_consistent(self, 'clean_bank')


class TestStep3MolID(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        PIPELINE.require('step3')
        with open(MOLID_INPUT) as f:
            cls.assignment = dict(line.rstrip('\n').rsplit(':', 1) for line in f if line.strip())
        with PIPELINE.log('step3_entities.log'):
            master = molidutils.pdb_to_masterlist(glob.glob(PIPELINE.path('clean_bank', '*.cif')))
        files_by_entity = molidutils.CreateMasterUniqueMolIDMapWithFileName(master)
        cls.files_by_entity = {e: {os.path.basename(f).split('+')[0] for f in files} for e, files in files_by_entity.items()}

    def test_outlier_has_entities_no_other_structure_has(self):
        only_outlier = [e for e, ids in self.files_by_entity.items() if ids == {OUTLIER}]
        self.assertTrue(only_outlier, '{0} no longer stands out as an outlier'.format(OUTLIER))

    def test_input_file_covers_entities_of_current_pdb_files(self):
        # Fails if the PDB renames an entity: the input file must then be updated.
        current = {e.strip() for e, ids in self.files_by_entity.items() if ids - {OUTLIER}}
        self.assertEqual(current, set(self.assignment))

    def test_outputs(self):
        subset = [a for a in ASSEMBLIES if not a.startswith(OUTLIER)]
        self.assertEqual(cif_files('standard_MolID_bank'), subset)
        csv = PIPELINE.path('standard_MolID_bank', 'OldChainID_NewChainID_numbered_MolID_Files_Map.csv')
        with open(csv) as f:
            rows = [line.strip().split(':') for line in f][1:]
        self.assertEqual({r[-1][:-4] for r in rows}, set(subset))

    def test_chain_ids_come_from_the_assignment(self):
        allowed = {c for chains in self.assignment.values() for c in chains.split(',')}
        for name in cif_files('standard_MolID_bank'):
            with self.subTest(structure=name):
                (asym,) = _support.cif_columns(PIPELINE.path('standard_MolID_bank', name + '.cif'),
                                               '_atom_site.label_asym_id')
                self.assertLessEqual(set(asym), allowed)

    def test_consistent_with_mmcif_dictionary(self):
        assert_consistent(self, 'standard_MolID_bank')


class TestStep4ChainID(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        PIPELINE.require('step4')
        cls.records = []
        with open(PIPELINE.path('standard_ChainID_bank', 'ChainStandardizationRecord.txt')) as f:
            for line in f:
                if line.strip():
                    path, old, new, score = line.strip().rsplit(':', 3)
                    cls.records.append((os.path.basename(path)[:-4], old, new, float(score)))

    def test_outputs(self):
        self.assertEqual(cif_files('standard_ChainID_bank'), cif_files('standard_MolID_bank'))
        self.assertEqual({r[0] for r in self.records}, set(cif_files('standard_MolID_bank')))

    def test_scores(self):
        for name, old, new, score in self.records:
            with self.subTest(structure=name, chain=old):
                self.assertGreaterEqual(score, 0.0)
                self.assertLessEqual(score, 1.0)
                if name.startswith('pdb_00001aw1'):
                    # 1aw1 has four assemblies, so it dominates the consensus.
                    self.assertEqual(score, 1.0)

    def test_different_organism_scores_lower(self):
        # 1btm (B. stearothermophilus) vs the 1aw1-dominated consensus: the outlier signal.
        best_1btm = max(r[3] for r in self.records if r[0].startswith('pdb_00001btm'))
        self.assertLess(best_1btm, 0.6)

    def test_renamed_chains_exist_in_output(self):
        for name, old, new, score in self.records:
            with self.subTest(structure=name, chain=old):
                (asym,) = _support.cif_columns(PIPELINE.path('standard_ChainID_bank', name + '.cif'),
                                               '_atom_site.label_asym_id')
                self.assertIn(new, asym)

    def test_standard_sequences_saved_as_json(self):
        # OCT2026: the consensus sequences are written as JSON, ready to be loaded with option 2) -> 5)
        with open(PIPELINE.path('standard_ChainID_bank', 'StandardSequences.json')) as f:
            std_seqs = json.load(f)
        # Every chain that Step 4 assigned is one of the standard sequences
        self.assertLessEqual({r[2] for r in self.records}, set(std_seqs))
        for chid, seq in std_seqs.items():
            with self.subTest(chain=chid):
                self.assertRegex(seq, r'^[A-Z]+$')

    def test_consistent_with_mmcif_dictionary(self):
        assert_consistent(self, 'standard_ChainID_bank')


class TestStep5ResidueID(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        PIPELINE.require('step5')
        cls.mapping = defaultdict(dict)  # structure -> {(chain, old number): (new number, ins code)}
        with open(PIPELINE.path('standard_ResidueID_bank', 'OldResID_NewResID_Map.csv')) as f:
            next(f)
            for line in f:
                old, new, cif = line.rstrip('\n').split(':')
                chain, oldnum = old.split('_', 1)
                parts = new.split()
                cls.mapping[cif[:-4]][(chain, oldnum.strip())] = (parts[0], parts[1] if len(parts) > 1 else '?')
        cls.structures = cif_files('standard_ResidueID_bank')
        cls.before = {s: residues(PIPELINE.path('standard_ChainID_bank', s + '.cif')) for s in cls.structures}
        cls.after = {s: residues(PIPELINE.path('standard_ResidueID_bank', s + '.cif')) for s in cls.structures}

    def test_outputs(self):
        self.assertEqual(self.structures, cif_files('standard_ChainID_bank'))
        self.assertEqual(set(self.mapping), set(self.structures))

    def test_cif_numbering_matches_map(self):
        # OCT2026: the new numbering is in auth_seq_id + pdbx_PDB_ins_code; label_seq_id is not changed
        for s in self.structures:
            with self.subTest(structure=s):
                self.assertEqual(len(self.before[s]), len(self.after[s]))
                for (chain, num, auth, ins, comp), (chain2, num2, auth2, ins2, comp2) in zip(self.before[s], self.after[s]):
                    self.assertEqual((chain, num, comp), (chain2, num2, comp2))
                    old = num + ('' if ins == '?' else ins)
                    self.assertEqual(self.mapping[s][(chain, old)], (auth2, ins2),
                                     'chain {0} residue {1}'.format(chain, old))

    def test_numbering_increases_along_each_chain(self):
        for s in self.structures:
            last = {}
            for chain, num, auth, ins, comp in self.after[s]:
                key = residue_order(auth, ins)
                if chain in last:
                    self.assertGreater(key, last[chain], '{0} chain {1} at {2}{3}'.format(s, chain, auth, ins))
                last[chain] = key

    def test_insertion_codes_assigned_for_gappy_columns(self):
        with_ins = [s for s in self.structures if any(r[3] != '?' for r in self.after[s])]
        self.assertTrue(with_ins, 'no insertion codes at occupancy threshold {0}'.format(OCCUPANCY_THRESHOLD))

    def test_same_number_same_residue_across_copies(self):
        # The four 1aw1 assemblies are copies of one protein: standardized numbering
        # must put the same residue type at the same number in all of them.
        by_position = defaultdict(set)
        for s in self.structures:
            if s.startswith('pdb_00001aw1'):
                for chain, num, auth, ins, comp in self.after[s]:
                    by_position[(chain, auth, ins)].add(comp)
        self.assertTrue(by_position)
        conflicts = {k: v for k, v in by_position.items() if len(v) > 1}
        self.assertFalse(conflicts)

    def test_consistent_with_mmcif_dictionary(self):
        assert_consistent(self, 'standard_ResidueID_bank')

    def test_no_alignment_files_in_banks(self):
        for bank in ('standard_ChainID_bank', 'standard_ResidueID_bank'):
            stray = [f for ext in ('fa', 'fasta', 'afa') for f in glob.glob(PIPELINE.path(bank, '*.' + ext))]
            self.assertFalse(stray)


if __name__ == '__main__':
    unittest.main()
