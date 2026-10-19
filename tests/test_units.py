"""
Unit tests for CIFcurator's pure functions. Offline, fast, no MUSCLE needed.

Run from the repository root:
    python -m unittest tests.test_units -v
"""
import contextlib
import io
import os
import pathlib
import shutil
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _support  # noqa: E402  (makes CIFcurator importable from src/)

from CIFcurator import alignmentutils, cifcurator_io, cleanutils, pdbutils  # noqa: E402
from CIFcurator import cifcuratorresiduestandardizationutils as resstd  # noqa: E402


def setUpModule():
    _support.quiet_known_warnings()


def quiet():
    """Silence the library's progress printing."""
    return contextlib.redirect_stdout(io.StringIO())


# Minimal mmCIF with every item cleanutils.simplify_cif reads.
# Two biological assemblies: assembly 1 = chain A (ALA, GLY), assembly 2 = chain B (SER, LYS, TRP).
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
ATOM 3 C CA . SER B 1 1 ? 1.000 9.000 1.000 1.00 10.00 1 SER B CA 1
ATOM 4 C CA . LYS B 1 2 ? 4.800 9.000 1.000 1.00 10.00 2 LYS B CA 1
ATOM 5 C CA . TRP B 1 3 ? 8.600 9.000 1.000 1.00 10.00 3 TRP B CA 1
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
_entity.id 1
_entity.type polymer
_entity.src_method man
_entity.pdbx_description 'Test protein'
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
    SEQRES = (">101m_A mol:protein length:154  MYOGLOBIN\n"
              "MVLSEGEWQLVLHVWAKVEAD\n"
              ">1aw1_A mol:protein length:255  TRIOSEPHOSPHATE ISOMERASE\n"
              "ARKFFVGGNWKMNGDKK\n"
              ">1aw1_B mol:protein length:255  TRIOSEPHOSPHATE ISOMERASE\n"
              "ARKFFVGGNWKMNGDKK\n"
              ">2y62_A mol:protein length:251  Triosephosphate isomerase\n"
              "MAPRKFFVGGNWKMNG\n")

    def test_get_idset_from_metadata(self):
        meta = ['>1aw1_B mol:protein', '>2y62_A mol:protein', '>1aw1_A mol:protein']
        self.assertEqual(pdbutils.get_idset_from_metadata(meta), ['1aw1', '2y62'])

    def test_keyword_search_is_case_insensitive_and_returns_next_line(self):
        seqfile = self.write('seqres.txt', self.SEQRES)
        seqs, meta = pdbutils.retrieve_sequence_from_PDB('triosephosphate isomerase', mode='metadata',
                                                         update=False, seqfile=seqfile)
        self.assertEqual(pdbutils.get_idset_from_metadata(meta), ['1aw1', '2y62'])
        self.assertEqual([s.strip() for s in seqs], ['ARKFFVGGNWKMNGDKK', 'ARKFFVGGNWKMNGDKK', 'MAPRKFFVGGNWKMNG'])

    def test_sequence_search_returns_preceding_header(self):
        seqfile = self.write('seqres.txt', self.SEQRES)
        seqs, meta = pdbutils.retrieve_sequence_from_PDB('FFVGGNWKMNG', mode='sequence',
                                                         update=False, seqfile=seqfile)
        self.assertEqual(len(seqs), 3)
        self.assertEqual(pdbutils.get_idset_from_metadata(meta), ['1aw1', '2y62'])

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

    def test_download_pdb_from_id_builds_rcsb_url(self):
        with mock.patch.object(pdbutils, 'download_from_url') as dl:
            pdbutils.download_pdb_from_id('1aw1', download_dir='/some/dir')
        dl.assert_called_once_with('https://files.rcsb.org/download/1AW1.cif', '/some/dir/1aw1.cif')


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


class TestSimplifyCif(TempDirTestCase):
    def simplify(self, text, bank):
        old = self.write('clean_bank/9xyz.cif', text)
        new = os.path.join(self.tmp, bank, '9xyz.cif')
        os.makedirs(os.path.dirname(new), exist_ok=True)
        cleanutils.simplify_cif(old, new, '.cif')
        return os.path.join(self.tmp, bank)

    def chains_in(self, path):
        (asym,) = _support.cif_columns(path, '_atom_site.label_asym_id')
        return sorted(set(asym))

    def test_one_file_per_assembly(self):
        bank = self.simplify(TWO_ASSEMBLIES, 'simple_bank')
        self.assertEqual(sorted(os.listdir(bank)), ['9xyz+01.cif', '9xyz+02.cif'])
        self.assertEqual(self.chains_in(os.path.join(bank, '9xyz+01.cif')), ['A'])
        self.assertEqual(self.chains_in(os.path.join(bank, '9xyz+02.cif')), ['B'])

    def test_single_assembly_gets_plus00(self):
        bank = self.simplify(ONE_ASSEMBLY, 'simple_bank')
        self.assertEqual(os.listdir(bank), ['9xyz+00.cif'])
        self.assertEqual(self.chains_in(os.path.join(bank, '9xyz+00.cif')), ['A', 'B'])

    def test_path_containing_cif_is_not_corrupted(self):
        # OCT2026 regression: ".cif" used to be applied as a regex to the whole path,
        # so "specific_cifbank" became "spec_bank" and the output landed elsewhere.
        bank = self.simplify(TWO_ASSEMBLIES, 'cifcurator_specific_cifbank')
        self.assertEqual(sorted(os.listdir(bank)), ['9xyz+01.cif', '9xyz+02.cif'])
        self.assertEqual(sorted(os.listdir(self.tmp)), ['cifcurator_specific_cifbank', 'clean_bank'])

    def test_output_is_readable_by_residue_step(self):
        bank = self.simplify(TWO_ASSEMBLIES, 'simple_bank')
        with quiet():
            seqs, resnums, structids, chids = resstd.pdb_to_structurelists(
                [os.path.join(bank, '9xyz+01.cif'), os.path.join(bank, '9xyz+02.cif')])
        self.assertEqual(sorted(seqs.values()), ['AG', 'SKW'])
        self.assertEqual(chids, ['A', 'B'])
        self.assertEqual(resnums[1]['B'], ['1 ', '2 ', '3 '])


class TestCleanCif(TempDirTestCase):
    def test_keeps_needed_categories_and_drops_others(self):
        text = TWO_ASSEMBLIES.replace('_exptl.method', '_struct_keywords.text FOO\n#\n_exptl.method')
        old = self.write('raw_bank/9xyz.cif', text)
        new = os.path.join(self.tmp, '9xyz_clean.cif')
        cleanutils.clean_cif(old, new)
        with open(new) as f:
            cleaned = f.read()
        self.assertTrue(cleaned.startswith('data_9XYZ'))
        self.assertNotIn('_struct_keywords', cleaned)
        for item in ('_entry.id', '_citation.title', '_entity.pdbx_description',
                     '_pdbx_struct_assembly_gen.asym_id_list', '_atom_site.group_PDB'):
            self.assertIn(item, cleaned)
        self.assertEqual(cleaned.count('\nATOM '), 5)


if __name__ == '__main__':
    unittest.main()
