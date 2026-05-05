'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

import os
import tempfile
import unittest

import numpy

import sasmol.config as config
import sasmol.system as system


class Test_unit_file_io_Files_biomt_metadata(unittest.TestCase):

    def _write_temp_pdb(self, lines):
        handle = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        try:
            handle.writelines(lines)
            return handle.name
        finally:
            handle.close()

    def test_read_pdb_without_biomt_records_empty_metadata(self):
        pdb_path = self._write_temp_pdb([
            "ATOM  53893  N   ILE N 515      73.944  41.799  41.652  1.00 36.37           N\n",
            "END\n",
        ])
        self.addCleanup(lambda: os.path.exists(pdb_path) and os.remove(pdb_path))

        molecule = system.Molecule(0)
        molecule.read_pdb(pdb_path)

        self.assertEqual(molecule.biomt(), {})

    def test_read_pdb_records_identity_biomt_without_coordinate_changes(self):
        pdb_path = self._write_temp_pdb([
            "REMARK 350 BIOMOLECULE: 1\n",
            "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC\n",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: N\n",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n",
            "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000\n",
            "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\n",
            "ATOM  53893  N   ILE N 515      73.944  41.799  41.652  1.00 36.37           N\n",
            "END\n",
        ])
        self.addCleanup(lambda: os.path.exists(pdb_path) and os.remove(pdb_path))

        molecule = system.Molecule(0)
        molecule.read_pdb(pdb_path)

        biomt = molecule.biomt()
        self.assertIn(1, biomt)
        self.assertEqual(biomt[1]['subdivs'], ['N'])
        self.assertEqual(biomt[1]['auth_bio_unit'], 'MONOMERIC')
        self.assertEqual(len(biomt[1]['rot']), 1)
        self.assertEqual(len(biomt[1]['trans']), 1)
        self.assertTrue(numpy.allclose(biomt[1]['rot'][0], numpy.identity(3)))
        self.assertTrue(numpy.allclose(biomt[1]['trans'][0], numpy.zeros(3)))

        expected = numpy.array([[[73.944, 41.799, 41.652]]], dtype=config.COORD_DTYPE)
        self.assertTrue(numpy.allclose(molecule.coor(), expected))
        self.assertTrue(any('REMARK 350   BIOMT1   1' in line for line in molecule.header()))

    def test_read_pdb_records_multiple_biomt_operators_and_chain_continuations(self):
        pdb_path = self._write_temp_pdb([
            "REMARK 350 BIOMOLECULE: 1\n",
            "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC\n",
            "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B\n",
            "REMARK 350                    AND CHAINS: C\n",
            "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000\n",
            "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000\n",
            "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000\n",
            "REMARK 350   BIOMT1   2  0.000000 -1.000000  0.000000       61.15700\n",
            "REMARK 350   BIOMT2   2 -1.000000  0.000000  0.000000       61.15700\n",
            "REMARK 350   BIOMT3   2  0.000000  0.000000 -1.000000       61.53850\n",
            "ATOM  53893  N   ILE A 515      73.944  41.799  41.652  1.00 36.37           N\n",
            "END\n",
        ])
        self.addCleanup(lambda: os.path.exists(pdb_path) and os.remove(pdb_path))

        molecule = system.Molecule(0)
        molecule.read_pdb(pdb_path)

        biomt = molecule.biomt()
        self.assertIn(1, biomt)
        self.assertEqual(biomt[1]['subdivs'], ['A', 'B', 'C'])
        self.assertEqual(biomt[1]['soft_bio_unit'], 'DIMERIC')
        self.assertEqual(len(biomt[1]['rot']), 2)
        self.assertEqual(len(biomt[1]['trans']), 2)
        self.assertAlmostEqual(float(biomt[1]['trans'][1][0]), 61.157, places=3)
        self.assertAlmostEqual(float(biomt[1]['trans'][1][2]), 61.5385, places=3)


if __name__ == '__main__':
    unittest.main()
