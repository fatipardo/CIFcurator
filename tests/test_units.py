"""
Unit tests for CIFcurator's pure functions. Offline, fast, no MUSCLE needed.

Run from the repository root:
    python -m unittest tests.test_units -v
"""
import contextlib
import gzip
import io
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _support  # noqa: E402  (makes CIFcurator importable from src/)

from CIFcurator import alignmentutils, cifcurator_io, cleanutils, pdbutils  # noqa: E402
from CIFcurator import cifcuratorresiduestandardizationutils as resstd  # noqa: E402
from CIFcurator import cifcuratorchainstandardizationutils as chainstd  # noqa: E402
from CIFcurator import cifcuratormolidcifutils as molidutils  # noqa: E402


def setUpModule():
    _support.quiet_known_warnings()


def quiet():
    """Silence the library's progress printing."""
    return contextlib.redirect_stdout(io.StringIO())


# Minimal mmCIF with every item cleanutils.clean_cif reads, and the categories the atoms refer to.
# Two biological assemblies: assembly 1 = chain A (ALA, GLY; SER not observed), assembly 2 = chain B (SER, LYS, TRP).
ATOM_HEADER = """loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""
ATOMS = """ATOM 1 C CA . ALA A 1 1 ? 1.000 1.000 1.000 1.00 10.00 1 ALA A CA 1
ATOM 2 C CA . GLY A 1 2 ? 4.800 1.000 1.000 1.00 10.00 2 GLY A CA 1
ATOM 3 C CA . SER B 2 1 ? 1.000 9.000 1.000 1.00 10.00 1 SER B CA 1
ATOM 4 C CA . LYS B 2 2 ? 4.800 9.000 1.000 1.00 10.00 2 LYS B CA 1
ATOM 5 N NZ . LYS B 2 2 ? 5.800 9.000 1.000 1.00 10.00 2 LYS B NZ 1
ATOM 6 C CA . TRP B 2 3 ? 8.600 9.000 1.000 1.00 10.00 3 TRP B CA 1
"""
COMMON = """data_9XYZ
#
_entry.id 9XYZ
#
_citation.title 'A synthetic test structure'
_citation.year 2026
_citation.pdbx_database_id_DOI ?
#
_citation_author.name 'Tester, T.'
#
_exptl.method 'X-RAY DIFFRACTION'
#
_refine.pdbx_refine_id 'X-RAY DIFFRACTION'
_refine.ls_d_res_high 2.00
#
loop_
_entity.id
_entity.type
_entity.src_method
_entity.pdbx_description
1 polymer man 'Test protein'
2 polymer man "Other: test protein's B chain"
#
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
2 polypeptide(L)
#
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 ALA n
1 2 GLY n
1 3 SER n
2 1 SER n
2 2 LYS n
2 3 TRP n
#
loop_
_chem_comp.id
_chem_comp.type
_chem_comp.mon_nstd_flag
_chem_comp.name
_chem_comp.formula
ALA 'L-peptide linking' y ALANINE 'C3 H7 N O2'
GLY 'peptide linking' y GLYCINE 'C2 H5 N O2'
LYS 'L-peptide linking' y LYSINE 'C6 H15 N2 O2 1'
SER 'L-peptide linking' y SERINE 'C3 H7 N O3'
TRP 'L-peptide linking' y TRYPTOPHAN 'C11 H12 N2 O2'
#
loop_
_struct_asym.id
_struct_asym.pdbx_blank_PDB_chainid_flag
_struct_asym.entity_id
A N 1
B N 2
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.ndb_seq_num
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.auth_seq_num
_pdbx_poly_seq_scheme.pdb_mon_id
_pdbx_poly_seq_scheme.auth_mon_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_ins_code
_pdbx_poly_seq_scheme.hetero
A 1 1 ALA 1 1 1 ALA ALA A . n
A 1 2 GLY 2 2 2 GLY GLY A . n
A 1 3 SER 3 3 ? ? ? A . n
B 2 1 SER 1 1 1 SER SER B . n
B 2 2 LYS 2 2 2 LYS LYS B . n
B 2 3 TRP 3 3 3 TRP TRP B . n
#
"""
TWO_ASSEMBLIES = COMMON + """loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A
2 1 B
#
""" + ATOM_HEADER + ATOMS + "#\n"
ONE_ASSEMBLY = COMMON + """_pdbx_struct_assembly_gen.assembly_id 1
_pdbx_struct_assembly_gen.oper_expression 1
_pdbx_struct_assembly_gen.asym_id_list A,B
#
""" + ATOM_HEADER + ATOMS + "#\n"
# A wwPDB assembly file, as downloaded in Step 1: no assembly tables, and the PDB ID is hidden.
ASSEMBLY_FILE = COMMON.replace('9XYZ', 'XXXX') + ATOM_HEADER + ATOMS + "#\n"


class TempDirTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='cifc_unit_')
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def write(self, relpath, text):
        path = os.path.join(self.tmp, relpath)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            f.write(text)
        return path


class TestResnConvert(unittest.TestCase):
    def test_standard_amino_acids(self):
        self.assertEqual(''.join(alignmentutils.ResnConvert(r) for r in ['MET', 'ALA', 'TRP', 'LYS']), 'MAWK')

    def test_modified_nucleotides_map_to_parent(self):
        self.assertEqual(alignmentutils.ResnConvert('PSU'), 'U')
        self.assertEqual(alignmentutils.ResnConvert('7MG'), 'G')

    def test_unknown_residue_is_X(self):
        self.assertEqual(alignmentutils.ResnConvert('HOH'), 'X')
        self.assertEqual(alignmentutils.ResnConvert('UNK'), 'X')


class TestGapPercentages(TempDirTestCase):
    def test_calculate_gap_percentages(self):
        gaps = alignmentutils.calculate_gap_percentages(['AC-D', 'A--D', 'ACGD', 'A--D'])
        self.assertEqual(list(gaps), [0.0, 50.0, 75.0, 0.0])

    def test_read_fasta_files(self):
        path = self.write('aln.fasta', '> Seq s1\nAC-\nDE\n> Seq s2\nA--DE\n')
        self.assertEqual(alignmentutils.read_fasta_files(path), ['AC-DE', 'A--DE'])


class TestScoreSequenceAlignment(unittest.TestCase):
    def test_identical(self):
        self.assertEqual(alignmentutils.ScoreSequenceAlignment('ACDE', 'ACDE'), 1.0)

    def test_fraction_of_reference_positions(self):
        self.assertEqual(alignmentutils.ScoreSequenceAlignment('AC-E', 'ACDF'), 0.5)

    def test_gaps_in_both_count_as_match(self):
        # Exact-identity scoring: '-' == '-' counts. Documented here so a change is noticed.
        self.assertEqual(alignmentutils.ScoreSequenceAlignment('A-', 'A-'), 1.0)


class TestPdbUtils(TempDirTestCase):
    # As in the wwPDB beta archive pdb_seqres.fasta.gz (extended PDB IDs)
    SEQRES = (">pdb_0000101m_A mol:protein length:154  MYOGLOBIN\n"
              "MVLSEGEWQLVLHVWAKVEAD\n"
              ">pdb_00001aw1_A mol:protein length:255  TRIOSEPHOSPHATE ISOMERASE\n"
              "ARKFFVGGNWKMNGDKK\n"
              ">pdb_00001aw1_B mol:protein length:255  TRIOSEPHOSPHATE ISOMERASE\n"
              "ARKFFVGGNWKMNGDKK\n"
              ">pdb_00002y62_A mol:protein length:251  Triosephosphate isomerase\n"
              "MAPRKFFVGGNWKMNG\n")
    IDS = ['pdb_00001aw1', 'pdb_00002y62']

    def test_extended_pdbid(self):
        for pdbid in ['1aw1', '1AW1', 'pdb_00001aw1', 'PDB_00001AW1', ' 1aw1\n']:
            with self.subTest(pdbid=pdbid):
                self.assertEqual(pdbutils.extended_pdbid(pdbid), 'pdb_00001aw1')
        for bad in ['1aw', 'pdb_1aw1', 'pdb_00001aw1x', '1aw1+01']:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                pdbutils.extended_pdbid(bad)

    def test_entry_url_hashes_on_characters_before_the_last(self):
        self.assertEqual(pdbutils.entry_url('1AW1'),
                         'https://files-beta.wwpdb.org/pub/wwpdb/pdb/data/entries/aw/pdb_00001aw1/')
        self.assertEqual(pdbutils.entry_url('pdb_1abc5678'),
                         'https://files-beta.wwpdb.org/pub/wwpdb/pdb/data/entries/67/pdb_1abc5678/')

    def test_get_idset_from_metadata(self):
        meta = ['>pdb_00001aw1_B mol:protein', '>pdb_00002y62_A mol:protein', '>pdb_00001aw1_A mol:protein']
        self.assertEqual(pdbutils.get_idset_from_metadata(meta), self.IDS)

    def test_get_idset_from_legacy_metadata(self):
        meta = ['>1aw1_B mol:protein', '>2y62_A mol:protein', '>pdb_00001aw1_A mol:protein']
        self.assertEqual(pdbutils.get_idset_from_metadata(meta), self.IDS)

    def test_keyword_search_is_case_insensitive_and_returns_next_line(self):
        seqfile = self.write('seqres.txt', self.SEQRES)
        seqs, meta = pdbutils.retrieve_sequence_from_PDB('triosephosphate isomerase', mode='metadata',
                                                         update=False, seqfile=seqfile)
        self.assertEqual(pdbutils.get_idset_from_metadata(meta), self.IDS)
        self.assertEqual([s.strip() for s in seqs], ['ARKFFVGGNWKMNGDKK', 'ARKFFVGGNWKMNGDKK', 'MAPRKFFVGGNWKMNG'])

    def test_sequence_search_returns_preceding_header(self):
        seqfile = self.write('seqres.txt', self.SEQRES)
        seqs, meta = pdbutils.retrieve_sequence_from_PDB('FFVGGNWKMNG', mode='sequence',
                                                         update=False, seqfile=seqfile)
        self.assertEqual(len(seqs), 3)
        self.assertEqual(pdbutils.get_idset_from_metadata(meta), self.IDS)

    def test_download_from_url(self):
        src = pathlib.Path(self.write('remote/file.cif', 'data_TEST\n'))
        target = os.path.join(self.tmp, 'local.cif')
        with quiet():
            pdbutils.download_from_url(src.as_uri(), target)
        with open(target) as f:
            self.assertEqual(f.read(), 'data_TEST\n')

    def test_download_from_url_rejects_truncated_response(self):
        class Short(io.BytesIO):
            def info(self):
                return {'Content-Length': '4096'}
        target = os.path.join(self.tmp, 'local.cif')
        with mock.patch.object(pdbutils, 'urlopen', return_value=Short(b'data_TEST\n')):
            with self.assertRaises(pdbutils.ContentTooShortError):
                pdbutils.download_from_url('https://example.org/file.cif', target)
        self.assertEqual(os.listdir(self.tmp), [])

    def test_download_from_url_concurrent_calls_do_not_share_staging_file(self):
        # While the first call is reading, a second, truncated download of the same target
        # runs and fails. The first must still publish its own complete data.
        target = os.path.join(self.tmp, 'local.cif')
        test = self

        class Response(io.BytesIO):
            def __init__(self, data, length):
                super().__init__(data)
                self.length = length

            def info(self):
                return {'Content-Length': str(self.length)}

        class Interrupted(Response):
            def read(self, *args):
                if self.tell() == 0:
                    with test.assertRaises(pdbutils.ContentTooShortError):
                        pdbutils.download_from_url('https://example.org/file.cif', target)
                return super().read(*args)

        full = b'data_TEST\n' * 100
        responses = iter([Interrupted(full, len(full)), Response(b'data_', len(full))])
        with mock.patch.object(pdbutils, 'urlopen', side_effect=lambda url: next(responses)), quiet():
            pdbutils.download_from_url('https://example.org/file.cif', target)
        with open(target, 'rb') as f:
            self.assertEqual(f.read(), full)
        self.assertEqual(os.listdir(self.tmp), ['local.cif'])

    def test_download_pdb_from_id_builds_beta_url(self):
        with mock.patch.object(pdbutils, 'download_from_url') as dl:
            pdbutils.download_pdb_from_id('1AW1', download_dir='/some/dir')
        dl.assert_called_once_with('https://files-beta.wwpdb.org/download/pdb_00001aw1.cif.gz',
                                   '/some/dir/pdb_00001aw1.cif', gunzip=True)

    def test_retrieve_seqfile_from_beta_archive(self):
        with mock.patch.object(pdbutils, 'download_from_url') as dl:
            pdbutils.retrieve_seqfile(seqfile='/some/dir/seqres.txt')
        dl.assert_called_once_with('https://files-beta.wwpdb.org/pub/wwpdb/pdb/derived_data/pdb_seqres.fasta.gz',
                                   '/some/dir/seqres.txt', gunzip=True)

    def test_download_from_url_gunzip(self):
        src = os.path.join(self.tmp, 'remote.cif.gz')
        with gzip.open(src, 'wb') as f:
            f.write(b'data_TEST\n')
        target = os.path.join(self.tmp, 'local.cif')
        with quiet():
            pdbutils.download_from_url(pathlib.Path(src).as_uri(), target, gunzip=True)
        with open(target) as f:
            self.assertEqual(f.read(), 'data_TEST\n')

    def test_list_assemblies_reads_entry_listing(self):
        listing = ('<a href="../">x</a><a href="pdb_00001aw1-assembly10.cif.gz">x</a>'
                   '<a href="pdb_00001aw1-assembly2.cif.gz">x</a><a href="pdb_00001aw1-assembly1.cif.gz">x</a>')
        with mock.patch.object(pdbutils, 'urlopen', return_value=io.BytesIO(listing.encode())) as op:
            names = pdbutils.list_assemblies_from_id('1AW1')
        op.assert_called_once_with('https://files-beta.wwpdb.org/pub/wwpdb/pdb/data/entries/aw/pdb_00001aw1/assemblies/')
        self.assertEqual(names, ['pdb_00001aw1-assembly1.cif.gz', 'pdb_00001aw1-assembly2.cif.gz',
                                 'pdb_00001aw1-assembly10.cif.gz'])

    def test_list_assemblies_of_entry_without_assemblies(self):
        missing = pdbutils.HTTPError('https://example.org/', 404, 'Not Found', {}, None)
        with mock.patch.object(pdbutils, 'urlopen', side_effect=missing):
            self.assertEqual(pdbutils.list_assemblies_from_id('pdb_00001aw1'), [])

    def test_download_assemblies_from_id_builds_beta_urls(self):
        names = ['pdb_00001aw1-assembly1.cif.gz', 'pdb_00001aw1-assembly2.cif.gz', 'pdb_00001aw1-assembly10.cif.gz']
        with mock.patch.object(pdbutils, 'list_assemblies_from_id', return_value=names), \
                mock.patch.object(pdbutils, 'download_from_url') as dl:
            pdbutils.download_assemblies_from_id('1aw1', download_dir='/some/dir')
        url = 'https://files-beta.wwpdb.org/pub/wwpdb/pdb/data/entries/aw/pdb_00001aw1/assemblies/'
        self.assertEqual(dl.call_args_list, [
            mock.call(url + 'pdb_00001aw1-assembly1.cif.gz', '/some/dir/pdb_00001aw1+01.cif', gunzip=True),
            mock.call(url + 'pdb_00001aw1-assembly2.cif.gz', '/some/dir/pdb_00001aw1+02.cif', gunzip=True),
            mock.call(url + 'pdb_00001aw1-assembly10.cif.gz', '/some/dir/pdb_00001aw1+10.cif', gunzip=True)])


class TestRunScript(TempDirTestCase):
    def test_interrupt_kills_script_and_its_children(self):
        # The script starts a long-running child, as the alignment scripts start MUSCLE.
        pidfile = os.path.join(self.tmp, 'child.pid')
        script = self.write('script.py', (
            'import subprocess, sys, time\n'
            'child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])\n'
            'open({0!r}, "w").write(str(child.pid))\n'
            'time.sleep(60)\n').format(pidfile))

        def interrupt(*args, **kwargs):
            for _ in range(100):
                if os.path.exists(pidfile) and os.path.getsize(pidfile):
                    break
                time.sleep(0.1)
            raise KeyboardInterrupt

        with mock.patch.object(subprocess.Popen, 'communicate', interrupt):
            with self.assertRaises(KeyboardInterrupt):
                _support.run_script(script, [], '', cwd=os.path.join(self.tmp, 'cwd'),
                                    log_path=os.path.join(self.tmp, 'log'))
        with open(pidfile) as f:
            child = f.read()
        # The orphaned child is reaped by init shortly after it is killed.
        for _ in range(50):
            state = subprocess.run(['ps', '-o', 'stat=', '-p', child],
                                   stdout=subprocess.PIPE, text=True).stdout.strip()
            if state in ('', 'Z'):
                break
            time.sleep(0.1)
        self.assertIn(state, ('', 'Z'), 'child process {0} still running'.format(child))


class TestCheckProject(TempDirTestCase):
    def test_create_and_delete_levels(self):
        proj = os.path.join(self.tmp, 'TIM')
        with quiet():
            cifcurator_io.check_project(projdir=proj)
            cifcurator_io.check_project(projdir=proj, level='raw_bank')
        self.assertTrue(os.path.isfile(os.path.join(proj, 'raw_bank', 'info.txt')))
        with quiet():
            cifcurator_io.check_project(projdir=proj, level='raw_bank', action='delete')
        self.assertFalse(os.path.exists(os.path.join(proj, 'raw_bank')))
        self.assertTrue(os.path.isdir(proj))


class TestCleanCif(TempDirTestCase):
    """Step 2 (clean, OCT2026: merged with simplify) on full PDB entries: one file per biological assembly."""
    def clean(self, text, bank):
        old = self.write('raw_bank/9xyz.cif', text)
        new = os.path.join(self.tmp, bank, '9xyz.cif')
        os.makedirs(os.path.dirname(new), exist_ok=True)
        cleanutils.clean_cif(old, new, '.cif')
        return os.path.join(self.tmp, bank)

    def chains_in(self, path):
        (asym,) = _support.cif_columns(path, '_atom_site.label_asym_id')
        return sorted(set(asym))

    def test_one_file_per_assembly(self):
        bank = self.clean(TWO_ASSEMBLIES, 'clean_bank')
        self.assertEqual(sorted(os.listdir(bank)), ['9xyz+01.cif', '9xyz+02.cif'])
        self.assertEqual(self.chains_in(os.path.join(bank, '9xyz+01.cif')), ['A'])
        self.assertEqual(self.chains_in(os.path.join(bank, '9xyz+02.cif')), ['B'])

    def test_single_assembly_gets_plus00(self):
        bank = self.clean(ONE_ASSEMBLY, 'clean_bank')
        self.assertEqual(os.listdir(bank), ['9xyz+00.cif'])
        self.assertEqual(self.chains_in(os.path.join(bank, '9xyz+00.cif')), ['A', 'B'])

    def test_path_containing_cif_is_not_corrupted(self):
        # OCT2026 regression: ".cif" used to be applied as a regex to the whole path,
        # so "specific_cifbank" became "spec_bank" and the output landed elsewhere.
        bank = self.clean(TWO_ASSEMBLIES, 'cifcurator_specific_cifbank')
        self.assertEqual(sorted(os.listdir(bank)), ['9xyz+01.cif', '9xyz+02.cif'])
        self.assertEqual(sorted(os.listdir(self.tmp)), ['cifcurator_specific_cifbank', 'raw_bank'])

    def test_output_is_readable_by_residue_step(self):
        bank = self.clean(TWO_ASSEMBLIES, 'clean_bank')
        with quiet():
            seqs, resnums, structids, chids = resstd.pdb_to_structurelists(
                [os.path.join(bank, '9xyz+01.cif'), os.path.join(bank, '9xyz+02.cif')])
        self.assertEqual(sorted(seqs.values()), ['AG', 'SKW'])
        self.assertEqual(chids, ['A', 'B'])
        self.assertEqual(resnums[1]['B'], ['1 ', '2 ', '3 '])

    def test_consistent_with_mmcif_dictionary(self):
        # OCT2026 reviewer: categories and key items required by the data in the file
        for text in (TWO_ASSEMBLIES, ONE_ASSEMBLY):
            bank = self.clean(text, 'clean_bank')
            for name in os.listdir(bank):
                with self.subTest(name=name):
                    self.assertEqual(_support.mmcif_problems(os.path.join(bank, name)), [])

    def test_keeps_only_this_assemblys_chains(self):
        bank = self.clean(TWO_ASSEMBLIES, 'clean_bank')
        asym, scheme, symbols = _support.cif_columns(os.path.join(bank, '9xyz+01.cif'), '_struct_asym.id',
                                                     '_pdbx_poly_seq_scheme.asym_id', '_atom_type.symbol')
        self.assertEqual(asym, ['A'])
        self.assertEqual(set(scheme), {'A'})
        self.assertEqual(symbols, ['C'])

    def test_drops_other_categories(self):
        text = TWO_ASSEMBLIES.replace('_exptl.method', '_struct_keywords.text FOO\n#\n_exptl.method')
        bank = self.clean(text, 'clean_bank')
        with open(os.path.join(bank, '9xyz+01.cif')) as f:
            cleaned = f.read()
        self.assertTrue(cleaned.startswith('data_9xyz+01\n'))
        self.assertNotIn('_struct_keywords', cleaned)
        self.assertNotIn('_pdbx_struct_assembly_gen', cleaned)
        self.assertEqual(cleaned.count('\nATOM '), 2)

    def test_entity_names_upper_case_without_quotes_or_colons(self):
        bank = self.clean(ONE_ASSEMBLY, 'clean_bank')
        (names,) = _support.cif_columns(os.path.join(bank, '9xyz+00.cif'), '_entity.pdbx_description')
        self.assertEqual(names, ['TEST PROTEIN', 'OTHER TEST PROTEINS B CHAIN'])


class TestCleanAssemblyFile(TempDirTestCase):
    """Step 2 on wwPDB assembly files (PDBID+NN.cif): each is rewritten whole, under the same name."""
    def clean_bank(self, names):
        for name in names:
            self.write('raw_bank/' + name, ASSEMBLY_FILE)
        os.makedirs(os.path.join(self.tmp, 'clean_bank'))
        with quiet():
            cleanutils.process(projdir=self.tmp, step='clean', source='raw_bank', target='clean_bank')
        return os.path.join(self.tmp, 'clean_bank')

    def test_one_output_per_assembly_file(self):
        bank = self.clean_bank(['9xyz+01.cif', '9xyz+02.cif', '8abc+01.cif'])
        self.assertEqual(sorted(os.listdir(bank)), ['8abc+01.cif', '9xyz+01.cif', '9xyz+02.cif'])

    def test_keeps_every_atom(self):
        bank = self.clean_bank(['9xyz+01.cif'])
        (atoms,) = _support.cif_columns(os.path.join(bank, '9xyz+01.cif'), '_atom_site.id')
        self.assertEqual(atoms, ['1', '2', '3', '4', '5', '6'])

    def test_entry_id_comes_from_filename(self):
        bank = self.clean_bank(['9xyz+01.cif'])
        with open(os.path.join(bank, '9xyz+01.cif')) as f:
            text = f.read()
        self.assertIn("_entry.id   '9XYZ'", text)
        self.assertNotIn('XXXX', text)
        (exptl,) = _support.cif_columns(os.path.join(bank, '9xyz+01.cif'), '_exptl.entry_id')
        self.assertEqual(exptl, ['9XYZ'])

    def test_extended_entry_id_stays_lower_case(self):
        # As in the wwPDB beta archive entries: _entry.id pdb_00001aw1
        bank = self.clean_bank(['pdb_00009xyz+01.cif'])
        with open(os.path.join(bank, 'pdb_00009xyz+01.cif')) as f:
            text = f.read()
        self.assertTrue(text.startswith('data_pdb_00009xyz+01\n'))
        self.assertIn("_entry.id   'pdb_00009xyz'", text)
        self.assertNotIn('XXXX', text)

    def test_consistent_with_mmcif_dictionary(self):
        bank = self.clean_bank(['9xyz+01.cif'])
        self.assertEqual(_support.mmcif_problems(os.path.join(bank, '9xyz+01.cif')), [])

    def test_simplify_is_the_same_step(self):
        self.write('raw_bank/9xyz+01.cif', ASSEMBLY_FILE)
        os.makedirs(os.path.join(self.tmp, 'clean_bank'))
        os.makedirs(os.path.join(self.tmp, 'simple_bank'))
        with quiet():
            cleanutils.process(projdir=self.tmp, step='clean', source='raw_bank', target='clean_bank')
            cleanutils.process(projdir=self.tmp, step='simplify', source='raw_bank', target='simple_bank')
        with open(os.path.join(self.tmp, 'clean_bank', '9xyz+01.cif')) as a, \
                open(os.path.join(self.tmp, 'simple_bank', '9xyz+01.cif')) as b:
            self.assertEqual(a.read(), b.read())


class TestStepsKeepFilesConsistent(TempDirTestCase):
    """
    OCT2026: Steps 3-5 rename chains and renumber residues. They must do the same in struct_asym and
    pdbx_poly_seq_scheme, so that the files stay consistent with the mmCIF dictionary.
    """
    def setUp(self):
        super().setUp()
        self.write('raw_bank/9xyz+01.cif', ASSEMBLY_FILE)
        os.makedirs(os.path.join(self.tmp, 'clean_bank'))
        os.makedirs(os.path.join(self.tmp, 'out'))
        with quiet():
            cleanutils.process(projdir=self.tmp, step='clean', source='raw_bank', target='clean_bank')
        self.path = os.path.join(self.tmp, 'clean_bank', '9xyz+01.cif')
        self.out = os.path.join(self.tmp, 'out', '9xyz+01.cif')

    def test_step3_renames_chains_everywhere(self):
        with quiet():
            master = molidutils.pdb_to_masterlist([self.path])
            master[0].chID_newchID_map.update({'A': 'X', 'B': 'Ab'})
            molidutils.masterlist_to_pdb([self.path], master, target_dir=os.path.join(self.tmp, 'out'))
        asym, auth, struct_asym, scheme = _support.cif_columns(
            self.out, '_atom_site.label_asym_id', '_atom_site.auth_asym_id', '_struct_asym.id',
            '_pdbx_poly_seq_scheme.pdb_strand_id')
        self.assertEqual(asym, ['X', 'X', 'Ab', 'Ab', 'Ab', 'Ab'])
        self.assertEqual(auth, asym)
        self.assertEqual(struct_asym, ['X', 'Ab'])
        self.assertEqual(scheme, ['X', 'X', 'X', 'Ab', 'Ab', 'Ab'])
        self.assertEqual(_support.mmcif_problems(self.out), [])

    def test_step3_concatenation_offsets_scheme_numbers(self):
        # Two entities concatenated into one chain C. mmCIF allows one entity per chain, so struct_asym keeps
        # a single row (the first entity's) and the file is not fully consistent: a known limitation.
        with quiet():
            master = molidutils.pdb_to_masterlist([self.path])
            master[0].chID_newchID_map.update({'A': 'C', 'B': 'C'})
            master[0].concat_order.update({'A': 1, 'B': 2})
            molidutils.masterlist_to_pdb([self.path], master, target_dir=os.path.join(self.tmp, 'out'))
        auth, struct_asym, strand, pdb_num = _support.cif_columns(
            self.out, '_atom_site.auth_seq_id', '_struct_asym.id', '_pdbx_poly_seq_scheme.pdb_strand_id',
            '_pdbx_poly_seq_scheme.pdb_seq_num')
        self.assertEqual(auth, ['1', '2', '50001', '50002', '50002', '50003'])
        self.assertEqual(struct_asym, ['C'])
        self.assertEqual(strand, ['C'] * 6)
        self.assertEqual(pdb_num, ['1', '2', '3', '50001', '50002', '50003'])
        self.assertEqual(_support.mmcif_problems(self.out),
                         ['_pdbx_poly_seq_scheme chain C entity 2 not in _struct_asym'] * 3 +
                         ['chain C: atoms in entity 2, _struct_asym says 1'] * 3)

    def test_step4_swaps_chains_everywhere(self):
        with quiet():
            chainstd.reassignedmaps_to_pdb([self.path], [{'A': 'B', 'B': 'A'}], [self.path],
                                           target_dir=os.path.join(self.tmp, 'out'))
        asym, struct_asym, scheme = _support.cif_columns(self.out, '_atom_site.label_asym_id', '_struct_asym.id',
                                                         '_pdbx_poly_seq_scheme.asym_id')
        self.assertEqual(asym, ['B', 'B', 'A', 'A', 'A', 'A'])
        self.assertEqual(struct_asym, ['B', 'A'])
        self.assertEqual(scheme, ['B', 'B', 'B', 'A', 'A', 'A'])
        self.assertEqual(_support.mmcif_problems(self.out), [])

    def test_step5_renumbers_author_numbering_only(self):
        template = {self.path: {'A_1 ': 10, 'A_2 ': '10 A', 'B_1 ': 1, 'B_2 ': 2, 'B_3 ': '2 A'}}
        resstd.conversiontemplate_to_pdb_FAPA([self.path], template, target_dir=os.path.join(self.tmp, 'out'))
        label, auth, ins = _support.cif_columns(self.out, '_atom_site.label_seq_id', '_atom_site.auth_seq_id',
                                                '_atom_site.pdbx_PDB_ins_code')
        self.assertEqual(label, ['1', '2', '1', '2', '2', '3'])
        self.assertEqual(auth, ['10', '10', '1', '2', '2', '2'])
        self.assertEqual(ins, ['?', 'A', '?', '?', '?', 'A'])
        pdb_num, auth_num, pdb_ins = _support.cif_columns(self.out, '_pdbx_poly_seq_scheme.pdb_seq_num',
                                                          '_pdbx_poly_seq_scheme.auth_seq_num',
                                                          '_pdbx_poly_seq_scheme.pdb_ins_code')
        # The unobserved SER (3rd row) keeps its number
        self.assertEqual(pdb_num, ['10', '10', '3', '1', '2', '2'])
        self.assertEqual(auth_num, ['10', '10', '?', '1', '2', '2'])
        self.assertEqual(pdb_ins, ['?', 'A', '.', '?', '?', 'A'])
        self.assertEqual(_support.mmcif_problems(self.out), [])


if __name__ == '__main__':
    unittest.main()
