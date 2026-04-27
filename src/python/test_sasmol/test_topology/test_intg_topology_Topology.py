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

import numpy
import sasmol.config as config
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

    def make_topology_test_molecule(self):
        molecule = system.Molecule(0)
        molecule.setNatoms(4)
        molecule.setNumber_of_frames(1)
        molecule.setAtom(['ATOM', 'ATOM', 'ATOM', 'HETATM'])
        molecule.setIndex(numpy.array([1, 2, 3, 4], numpy.int32))
        molecule.setName(['N', 'H', 'P', 'O'])
        molecule.setLoc([' ', ' ', ' ', ' '])
        molecule.setResname(['ALA', 'ALA', 'ADE', 'HOH'])
        molecule.setResid(numpy.array([1, 1, 2, 3], numpy.int32))
        molecule.setChain(['A', 'A', 'B', 'W'])
        molecule.setCoor(numpy.zeros((1, 4, 3), config.COORD_DTYPE))
        molecule.setRescode([' ', ' ', ' ', ' '])
        molecule.setSegname(['SEG1', 'SEG1', 'SEG2', 'WAT'])
        molecule.setElement(['N', 'H', 'P', 'O'])
        molecule.setCharge([' ', ' ', ' ', ' '])
        molecule.setMoltype(['protein', 'protein', 'nucleic', 'water'])
        molecule.setOccupancy(['0.25', '0.25', '0.25', '0.25'])
        molecule.setBeta(['0.50', '0.50', '0.50', '0.50'])
        return molecule

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

    def test_make_constraint_pdb_sets_backbone_occupancy(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '1CRN.pdb')
        outfile = self.make_temp_pdb()

        molecule.make_constraint_pdb(outfile, 'backbone', field='occupancy')

        constrained = system.Molecule(0)
        constrained.read_pdb(outfile)

        self.assertEqual(constrained.occupancy()[0], '1.00')
        self.assertEqual(constrained.occupancy()[1], '1.00')
        self.assertEqual(constrained.occupancy()[4], '0.00')

    def test_make_constraint_pdb_sets_heavy_atoms(self):
        molecule = self.make_topology_test_molecule()
        outfile = self.make_temp_pdb()

        molecule.make_constraint_pdb(outfile, 'heavy')

        constrained = system.Molecule(0)
        constrained.read_pdb(outfile)

        self.assertEqual(constrained.beta()[0], '1.00')
        self.assertEqual(constrained.beta()[1], '0.00')
        self.assertEqual(constrained.beta()[2], '1.00')

    def test_make_constraint_pdb_sets_protein_atoms(self):
        molecule = self.make_topology_test_molecule()
        outfile = self.make_temp_pdb()

        molecule.make_constraint_pdb(outfile, 'protein')

        constrained = system.Molecule(0)
        constrained.read_pdb(outfile)

        self.assertEqual(constrained.beta()[0], '1.00')
        self.assertEqual(constrained.beta()[1], '1.00')
        self.assertEqual(constrained.beta()[2], '0.00')
        self.assertEqual(constrained.beta()[3], '0.00')

    def test_make_constraint_pdb_sets_nucleic_atoms(self):
        molecule = self.make_topology_test_molecule()
        outfile = self.make_temp_pdb()

        molecule.make_constraint_pdb(outfile, 'nucleic')

        constrained = system.Molecule(0)
        constrained.read_pdb(outfile)

        self.assertEqual(constrained.beta()[0], '0.00')
        self.assertEqual(constrained.beta()[2], '1.00')
        self.assertEqual(constrained.beta()[3], '0.00')

    def test_make_constraint_pdb_sets_solute_atoms(self):
        molecule = self.make_topology_test_molecule()
        outfile = self.make_temp_pdb()

        molecule.make_constraint_pdb(outfile, 'solute')

        constrained = system.Molecule(0)
        constrained.read_pdb(outfile)

        self.assertEqual(constrained.beta()[0], '1.00')
        self.assertEqual(constrained.beta()[1], '1.00')
        self.assertEqual(constrained.beta()[2], '1.00')
        self.assertEqual(constrained.beta()[3], '0.00')

    def test_make_constraint_pdb_reset_false_preserves_unselected_values(self):
        molecule = self.make_topology_test_molecule()
        outfile = self.make_temp_pdb()

        molecule.make_constraint_pdb(outfile, 'nucleic', reset=False)

        constrained = system.Molecule(0)
        constrained.read_pdb(outfile)

        self.assertEqual(constrained.beta()[0], '0.50')
        self.assertEqual(constrained.beta()[2], '1.00')
        self.assertEqual(constrained.beta()[3], '0.50')

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

    def test_make_backbone_pdb_from_fasta_uses_n_terminal_glycine_patch(self):
        molecule = system.Molecule(0)
        molecule.setFasta(['G', 'A'])
        outfile = self.make_temp_pdb()

        molecule.make_backbone_pdb_from_fasta(outfile, 'protein')

        backbone = system.Molecule(0)
        backbone.read_pdb(outfile)

        self.assertEqual(backbone.resname()[0], 'GLYP')
        self.assertEqual(backbone.resname()[1], 'ALA')

    def test_make_backbone_pdb_from_fasta_uses_n_terminal_proline_patch(self):
        molecule = system.Molecule(0)
        molecule.setFasta(['P', 'A'])
        outfile = self.make_temp_pdb()

        molecule.make_backbone_pdb_from_fasta(outfile, 'protein')

        backbone = system.Molecule(0)
        backbone.read_pdb(outfile)

        self.assertEqual(backbone.resname()[0], 'PROP')
        self.assertEqual(backbone.resname()[1], 'ALA')

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

    def test_create_fasta_format_by_segname(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '3AAD-2chain.pdb')
        molecule.create_fasta(fasta_format=True, by_segname=True, name='demo')

        result = molecule.fasta()
        self.assertIn('>demo segname:M', result)
        self.assertIn('>demo segname:N', result)

    def test_create_fasta_format_width(self):
        molecule = system.Molecule(0)
        molecule.read_pdb(DataPath + '1CRN.pdb')
        molecule.create_fasta(fasta_format=True, width='10', name='demo')

        result = molecule.fasta()
        self.assertIn('TTCCPSIVAR\nSNFNVCRLPG', result)

    def test_create_fasta_excludes_hetatm_when_requested(self):
        molecule = self.make_topology_test_molecule()

        molecule.create_fasta()
        self.assertEqual(molecule.fasta(), ['A', 'A', 'X'])

        molecule.create_fasta(fasta_format=True, exclude_hetatm=True)
        self.assertEqual(molecule.fasta(), '>\nAA\n')

    def tearDown(self):
        for filename in self.tempfiles:
            if os.path.exists(filename):
                os.remove(filename)


if __name__ == '__main__':
    main()
