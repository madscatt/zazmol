#import unittest
#from operate import Move  # Adjust the import based on your project structure
#from sasmol import system  # Assuming sasmol.system.Molecule is used in Move



from sasmol.test_sasmol.utilities import env,util
from unittest import main, skipIf
import unittest
import sasmol.system as system
import sasmol.operate as operate

import numpy

import warnings

import os
floattype=os.environ['SASMOL_FLOATTYPE']



DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep


class TestMove(unittest.TestCase):

    def setUp(self):
        # Initialize the Move object and any other necessary objects
#        self.move = Move()
        self.molecule_1 = system.Molecule()
        self.molecule_2 = system.Molecule()
        
        # Load or create test data for the molecules
        # For example, you might load a PDB file or create a mock molecule
        self.molecule_1.read_pdb(DataPath+'1CRN.pdb')
        self.molecule_2.read_pdb(DataPath+'1CRN.pdb')

        axis = 'z'
        frame = 0
        theta=numpy.pi/2.0
        #
        self.molecule_2.rotate(frame,axis,theta)

    def test_align_trajectory(self):
        # Define the basis and frame for the alignment
        self_basis = 'CA'
        other_basis = 'CA'
        frame = 0

        # Initialize the subsets
        subset_self, subset_other = self.molecule_1.initialize_subsets(self.molecule_1, self_basis, other_basis, frame)

        # Call the align_trajectory method with the precomputed subsets
        self.align_trajectory(self.molecule_1, self.molecule_2, self_basis, other_basis, subset_self=subset_self, subset_other=subset_other, frame=frame)

        # Add assertions to verify the expected behavior
        # For example, you might check the coordinates of the aligned molecules
        aligned_coords_1 = self.molecule_1.coor()[frame]
        aligned_coords_2 = self.molecule_2.coor()[frame]

        # Example assertion (adjust based on your specific requirements)
        self.assertTrue((aligned_coords_1 == aligned_coords_2).all(), "The molecules are not aligned correctly")

if __name__ == '__main__':
    unittest.main()
