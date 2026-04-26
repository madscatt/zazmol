'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from sasmol.test_sasmol.utilities import env

from unittest import main
import tempfile
import unittest
import os

import sasmol.system as system


DataPath = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    '..',
    'data',
    'pdb_common'
) + os.path.sep


class Test_intg_topology_Topology(unittest.TestCase):

    def setUp(self):
        self.tempfiles = []

    def make_temp_pdb(self):
        handle = tempfile.NamedTemporaryFile(
            suffix='.pdb',
            delete=False,
            dir=os.getcwd()
        )
        handle.close()
        self.tempfiles.append(handle.name)
        return handle.name

    def test_renumber_default(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '2AAD.pdb')
        molecule.renumber(index=10, resid=20)
        molecule.renumber()

        self.assertEqual(molecule.index()[0], 1)
        self.assertEqual(molecule.index()[-1], molecule.natoms())
        self.assertEqual(molecule.resid()[0], 1)
        self.assertEqual(molecule.resid()[-1], 2)

    def test_renumber_index_only(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '1ATM.pdb')
        expected_resid = molecule.resid()[0]
        molecule.renumber(index=3)

        self.assertEqual(molecule.index()[0], 3)
        self.assertEqual(molecule.resid()[0], expected_resid)

    def test_renumber_resid_only(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '2AAD.pdb')
        molecule.renumber(resid=8)

        self.assertEqual(molecule.index()[0], 1)
        self.assertEqual(molecule.resid()[0], 8)
        self.assertEqual(molecule.resid()[-1], 9)

    def test_make_constraint_pdb_sets_backbone_beta(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '1CRN.pdb')
        outfile = self.make_temp_pdb()

        molecule.make_constraint_pdb(outfile, 'backbone')

        constrained = system.Molecule(0)
        constrained.read_pdb(outfile)

        self.assertEqual(constrained.beta()[0], '1.00')
        self.assertEqual(constrained.beta()[1], '1.00')
        self.assertEqual(constrained.beta()[4], '0.00')

    def test_make_backbone_pdb_from_fasta_protein(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '1CRN.pdb')
        molecule.create_fasta()
        outfile = self.make_temp_pdb()

        molecule.make_backbone_pdb_from_fasta(outfile, 'protein')

        backbone = system.Molecule(0)
        backbone.read_pdb(outfile)

        self.assertEqual(backbone.natoms(), len(molecule.fasta()))
        self.assertEqual(backbone.name()[0], 'CA')

    def test_make_backbone_pdb_from_fasta_nucleic(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + 'rna.pdb')
        molecule.create_fasta()
        outfile = self.make_temp_pdb()

        molecule.make_backbone_pdb_from_fasta(outfile, 'nucleic')

        backbone = system.Molecule(0)
        backbone.read_pdb(outfile)

        self.assertEqual(backbone.natoms(), len(molecule.fasta()))
        self.assertEqual(backbone.name()[0], "O5'")

    def test_create_fasta_default(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '1CRN.pdb')
        molecule.create_fasta()

        self.assertEqual(molecule.fasta()[:5], ['T', 'T', 'C', 'C', 'P'])

    def test_create_fasta_format_by_chain(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '3AAD-2chain.pdb')
        molecule.create_fasta(fasta_format=True, by_chain=True, name='demo')

        result = molecule.fasta()
        self.assertIn('>demo chain:M', result)
        self.assertIn('>demo chain:N', result)

    def tearDown(self):
        for filename in self.tempfiles:
            if os.path.exists(filename):
                os.remove(filename)


if __name__ == '__main__':
    main()
