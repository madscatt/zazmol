from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
#
#
#    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	SUBSET
#
#	01/04/2011	--	initial coding 			            :	jc
#	08/19/2016	--	added doc strings                   :	jc
#
# LC	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
import os
import sys
import string
import locale
import struct
import numpy
import time
import sasmol.mask
import sasmol.util as utilities
import random

'''
	Subset is the module that contains the classes that 
	allow users to extract a subset of objects or values from 
	objects.  These subsets are used to do things like align
	molecules, check for overlap, filter molecular data to
	select for fields in the PDB file.	

	These classes are accessed by the Atom class found in
	the system module.

'''


class Mask(object):

    """ Base class containing methods to extract or combine system objects
        using numpy masks

        Examples
        ========

        First example shows how to use class methods from system object:

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> basis_filter = 'name[i] == "CA" and resid[i] < 10'
        >>> error, mask = molecule.get_subset_mask(basis_filter)  
        >>> import numpy
        >>> numpy.nonzero(mask)
        (array([  4,  11,  21,  45,  55,  66,  82, 101, 112]),) 

        Note
        ----

        `self` parameter is not shown in the ``Parameters`` section in the documentation

    """

    def get_dihedral_subset_mask(self, flexible_residues, mtype):
        '''
        This method creates an array of ones and/or zeros
        of the length of the number of atoms in "self". It
        uses the user-supplied flexible_residue array to
        determine which atoms to include in the mask.
        This version is hard-wired for proteins or rna to choose
        the C(n-1), N(n), CA(n), C(n), and N(n+1) atoms		
        or the O3'(n-1), P(n), O5'(n), C5'(n), C4'(n), 
        C3'(n), O3'(n) and P(n+1) atoms that form the basis set		
        for the rotation phi & psi or alpha, beta, delta, epsilon, and eta
        angles respectively.  This method calles a c-method called
        mask to speed up the calculation (24.5 X faster).	
        '''
        natoms = self.natoms()
        name = self.name()
        resid = self.resid().astype(numpy.long)
        nflexible = len(flexible_residues)

        nresidues = int(self.resid()[-1] - self.resid()[0] + 1)

        farray = numpy.zeros((nflexible, natoms), numpy.long)

        sasmol.mask.get_mask_array(farray, name, resid,
                                   flexible_residues, nresidues, mtype)

        return farray

    def init_child(self, descriptor):
        '''
        This method allows one to create a list of Molecule objects
        that are defined by the input descriptor.


        usage:

                This is a way to create a mask to be used somewhere else:

                m1=system.Molecule(0)	### create a molecule m1
                m1.read_pdb(filename)	### read in variables, coor, etc.
                m1.initialize_children()   ### set up the masks etc.

                . . . do stuff . . . 

                This initializes the following "children" with their
                masks already defined to the "parent" molecule

                names() 		: names_mask()
                resnames()		: resnames_mask()
                resids()		: resids_mask()
                chains()		: chains_mask()
                segnames()		: segnames_mask()
                occupancies()		: occupancies_mask()
                betas()			: betas_mask()
                elements()		: elements_mask()

                The objects on the left contain the unique values and
                the objects on the right contain the masks that have
                the indices to extract the information for each unique
                value from the parent molecule.

                NOTE: the pluarity of the words is chosen for a reason
                to distinguish the singular words used to keep track
                of the parent variables (name --> name[i] for each
                atom, while names --> corresponds to the unique names
                in the parent: len(names) <= len(name))

                For "min3.pdb" if one wants to know the unique elements
                you would type:

                m1.elements()

                which yields:

                ['N', 'H', 'C', 'O', 'S', 'ZN']

                So, given a pre-defined object that has atomic information
                initialized by reading in the PDB file and intializing all
                children as shown above, one can get a list of subset
                objects for each type of element by typing:

                element_molecules = m1.init_child('elements')

                then you could parse the full-subset molecule as its own
                entity

                com = element_molecules[0].calccom(0)

                which would give the center of mass for all the "N" atoms
                in the parent molecule.

                Another example would be to get the COM of each amino acid
                in a protein.

                residue_molecules = m1.init_child('resids')

                for i in xrange(m1.number_of_resids()):
                        print(residue_molecules[i].calccom(0))

                NOTE: coordinates will have to be updated separately using
                        get_coor_using_mask ... using the mask(s) generated
                        by file_io.initialize_children()

        '''

        import sasmol.system as system

        at = 'len(self.' + descriptor + '())'
        number_of_objects = eval(at)

        frame = 0
        object_list = []

        for i in xrange(number_of_objects):

            new_object = system.Molecule(0)
            at = 'self.' + descriptor + '_mask()[' + str(i) + ']'
            mask = eval(at)
            error = self.copy_molecule_using_mask(new_object, mask, frame)

            object_list.append(new_object)

        return object_list

    def get_subset_mask(self, basis_filter):
        '''
        This method creates an array of ones and/or zeros
        of the length of the number of atoms in "self" and
        uses the user-supplied filter string to filter the
        parameter descriptors to obtain a subset array
        that can be used to filter entities in other methods
        either in this class or elsewhere.

        usage:

                Here is a way to create a mask to be used somewhere else:

                m1=system.Molecule(0)	### create a molecule m1
                m1.read_pdb(filename)	### read in variables, coor, etc.

                . . . do stuff . . . 

                basis_filter = XXXX     ### your (see examples below)

                error,mask = m1.get_subset_mask(basis_filter)  ### get a mask

                . . . do something with the mask using other functions in this class . . . 

                Here are some example basis_filter strings:

                basis_filter = 'name[i] == "CA" and resid[i] < 10'
                basis_filter = 'name[i][0] == "H" and resid[i] < 10'
                basis_filter = 'name[i] == "CA" and resid[i] >= 1 and resid[i] < 10'

                The syntax for basis selection can be quite eloborate.  For example,

                basis_filter = 'name[i] == "CA" and resid[i] >= 1 and resid[i] < 10 and moltype=="protein" and chain=="F" and occupancy==1 and beta>10.0 and element=="C" ...'

                could be used for advanced selection needs.  See API for full details.

        '''
        index = self.index()
        name = self.name()
        loc = self.loc()
        resname = self.resname()
        chain = self.chain()
        resid = self.resid()
        rescode = self.rescode()
        occupancy = self.occupancy()
        beta = self.beta()
        segname = self.segname()
        element = self.element()
        charge = self.charge()
        moltype = self.moltype()
        residue_flag = self.residue_flag()

#
# OPEN	Need to add other properties (surface, buried, helix, sheet, turn, polar
#
# OPEN	Need to add other properties (hydrophobic, etc.)
#
# OPEN	Need to filter string and replace certain keywords with atom names (calpha, backbone)
#
# OPEN	Need to handle "all" keyword as well
#
# OPEN	Need to add a distance selection logic: find atoms greater than 5 angstroms from ...
#
#
#	basis_filter == 'name[i] == "CA" less_than 5 angstroms from name[j] == "N" '
#
#	basis_filter == 'name[i] == "CA" greater_than 5 angstroms from name[j] == "N" '
#
#	basis_filter == 'name[i] == "CA" between 3 and 5 angstroms from name[j] == "N" '
#
#	ANGLES
#
#	basis_filter == '(segname[i] == "WAT" and (name[i] == "OH2" or name[i] == "H1")) less_than 104 degrees from (segname[j] == "WAT" and name[j] == "H1") '
#
#	basis_filter == '(segname[i] == "WAT" and name[i] == "OH2") between 100 and 106 degrees from (segname[j] == "WAT" and name[j] == "H1") '
#
# 	--- or have the code pull out the donors and acceptors from selections?
#
#	basis_filter == '(resid[i] == "TIP3") between 100 and 106 degrees from (resid[j] == "TIP3")'
#
#	--- or define a hydrogen_bonded keyword that takes care of distance/angle & donor/acceptor
#
#	basis_filter == '(resid[i] == "TIP3") hydrogen_bonded to (resid[j] == "TIP3")'
#	basis_filter == '(segment[i] == "RNA") hydrogen_bonded to (resid[j] == "TIP3")'
#	basis_filter == '(segment[i] == "RNA") hydrogen_bonded to (segment[j] == "RNA")'
#
#	--- or define an atom as a donor or acceptor and use hydrogen_bonded keyword for distance/angle
#
#	basis_filter == '(segment[i] == "RNA" and donor[i] == 1) hydrogen_bonded to (resid[j] == "TIP3")'
#	basis_filter == '(segment[i] == "RNA" and acceptor[i] == 1) hydrogen_bonded to (resid[j] == "TIP3")'
#
#
        error = []

        preliminary_mask_array = []

        natoms = self.natoms()

        for i in xrange(natoms):

            try:
                if(eval(basis_filter)):
                    preliminary_mask_array.append(1)
                else:
                    preliminary_mask_array.append(0)
            except:
                error.append('failed to evaluate filter selection ' +
                             basis_filter + ' for atom ' + str(i))
                return error, preliminary_mask_array

        if(numpy.sum(preliminary_mask_array) == 0):
            error.append(
                'found no atoms using filter selection ' + basis_filter)
            return error, preliminary_mask_array
        else:
            mask_array = numpy.array(preliminary_mask_array, numpy.int32)

        return error, mask_array
        '''


		natoms = self.natoms()
		mask_array = numpy.zeros(natoms,numpy.int32)

		for i in xrange(natoms):
			try:
				if(eval(basis_filter)):
					mask_array[i]=1
			except:
				error.append('failed to evaluate filter selection '+basis_filter+' for atom '+str(i))
				return error, mask_array		
		if(numpy.sum(mask_array) == 0):
				error.append('found no atoms using filter selection '+basis_filter)
				return error, mask_array		
		#if not qfound:
		#		error.append('found no atoms using filter selection '+basis_filter)
		#		return error, mask_array	
		return error, mask_array
		'''

    def merge_two_molecules(self, mol1, mol2):
        '''
        This method combines two molecules into a single, new molecule. 
        It will assign coordinates from the first frame of a molecule.

        usage:


                m1=system.Molecule(0)	### create a molecule m1
                m1.read_pdb(filename1)	### read in variables, coor, etc.

                m2=system.Molecule(1)	### create a molecule m2
                m2.read_pdb(filename2)	### read in variables, coor, etc.

                m3=system.Molecule(2)	### create a molecule m3


                . . . do stuff . . . 

                error = m3.merge_two_molecules(m1,m2) 	### sets the values that define mol3

        '''
        error = []
        atom = []
        index = []
        name = []
        loc = []
        resname = []
        chain = []
        resid = []
        rescode = []
        x = []
        y = []
        z = []
        occupancy = []
        beta = []
        segname = []
        element = []
        charge = []
        moltype = []
        original_index = []
        original_resid = []
        residue_flag = []

        natoms1 = mol1.natoms()
        natoms2 = mol2.natoms()

        print('natoms1 = ', natoms1)
        print('natoms2 = ', natoms2)

        frame = 0

        for i in xrange(natoms1):
            try:
                # if True:
                atom.append(mol1._atom[i])
                index.append(mol1._index[i])
                name.append(mol1._name[i])
                loc.append(mol1._loc[i])
                resname.append(mol1._resname[i])
                chain.append(mol1._chain[i])
                resid.append(mol1._resid[i])
                rescode.append(mol1._rescode[i])
                x.append(mol1.coor()[frame, i, 0])
                y.append(mol1.coor()[frame, i, 1])
                z.append(mol1.coor()[frame, i, 2])
                occupancy.append(mol1._occupancy[i])
                beta.append(mol1._beta[i])
                segname.append(mol1._segname[i])
                element.append(mol1._element[i])
                charge.append(mol1._charge[i])
                moltype.append(mol1._moltype[i])
                original_index.append(mol1._original_index[i])
                original_resid.append(mol1._original_resid[i])
                residue_flag.append(mol1._residue_flag[i])

            except:
                error.append(
                    'failed in copy_molecule when attempting to assign descriptors to atom ' + str(i) + ' from mol1')
                return error

        last_index_mol1 = mol1._index[-1]
        this_index = last_index_mol1 + 1

        for i in xrange(natoms2):
            try:
                # if True:
                atom.append(mol2._atom[i])
                index.append(this_index)
                name.append(mol2._name[i])
                loc.append(mol2._loc[i])
                resname.append(mol2._resname[i])
                chain.append(mol2._chain[i])
                resid.append(mol2._resid[i])
                rescode.append(mol2._rescode[i])
                x.append(mol2.coor()[frame, i, 0])
                y.append(mol2.coor()[frame, i, 1])
                z.append(mol2.coor()[frame, i, 2])
                occupancy.append(mol2._occupancy[i])
                beta.append(mol2._beta[i])
                segname.append(mol2._segname[i])
                element.append(mol2._element[i])
                charge.append(mol2._charge[i])
                moltype.append(mol2._moltype[i])
                original_index.append(mol2._original_index[i])
                original_resid.append(mol2._original_resid[i])
                residue_flag.append(mol2._residue_flag[i])
                this_index = this_index + 1

            except:
                error.append(
                    'failed in copy_molecule when attempting to assign descriptors to atom ' + str(i) + ' from mol2')
                return error
        
        x = numpy.array(x, numpy.float)
        y = numpy.array(y, numpy.float)
        z = numpy.array(z, numpy.float)

        coor = numpy.zeros((1, natoms1 + natoms2, 3), numpy.float)

        try:
            coor[frame, :, 0] = x
            coor[frame, :, 1] = y
            coor[frame, :, 2] = z
        except:
            error.append('failed in merge molecule when assigning coordinates')
            return error

        self._atom = atom
        self._index = index
        self._name = name
        self._loc = loc
        self._resname = resname
        self._chain = chain
        self._resid = resid
        self._rescode = rescode
        self._occupancy = occupancy
        self._beta = beta
        self._segname = segname
        self._element = element
        self._charge = charge
        self._moltype = moltype
        self._coor = numpy.array(coor)
        self._natoms = len(index)
        self._original_index = original_index
        self._original_resid = original_resid
        self._residue_flag = residue_flag

        self._unique_names = list(numpy.unique(self._name))
        self._unique_resnames = list(numpy.unique(self._resname))
        self._unique_resids = list(numpy.unique(self._resid))
        self._unique_chains = list(numpy.unique(self._chain))
        self._unique_segnames = list(numpy.unique(self._segname))
        self._unique_occupancies = list(numpy.unique(self._occupancy))
        self._unique_betas = list(numpy.unique(self._beta))
        self._unique_elements = list(numpy.unique(self._element))
        self._unique_moltypes = list(numpy.unique(self._moltype))

        self._number_of_names = len(self._unique_names)
        self._number_of_resnames = len(self._unique_resnames)
        self._number_of_resids = len(self._unique_resids)
        self._number_of_chains = len(self._unique_chains)
        self._number_of_segnames = len(self._unique_segnames)
        self._number_of_occupancies = len(self._unique_occupancies)
        self._number_of_betas = len(self._unique_betas)
        self._number_of_elements = len(self._unique_elements)
        self._number_of_moltypes = len(self._unique_moltypes)

        self._conect = mol1._conect

        try:
            for ndx, list_ndxs in mol2._conect.iteritems():
                self._conect[ndx] = list_ndxs
        except:
            pass

        return error

    def copy_molecule_using_mask(self, other, mask, frame):
        '''
        This method initializes the standard descriptors and
        coordinates for a subset molecule defined by the
        supplied mask array.

        usage:

                Here is a way to create a mask to be used somewhere else:

                m1=system.Molecule(0)	### create a molecule m1
                m1.read_pdb(filename)	### read in variables, coor, etc.

                . . . do stuff . . . 

                basis_filter = XXXX     ### your (see examples below)

                error,mask = m1.get_subset_mask(basis_filter)  ### get a mask

                sub_m1=system.Molecule(1)		### create a new molecule sub_m1

                error = m1.copy_molecule_using_mask(sub_m1,mask,frame) ### initializes sub_m1

        '''
        error = []
        atom = []
        index = []
        name = []
        loc = []
        resname = []
        chain = []
        resid = []
        rescode = []
        occupancy = []
        beta = []
        segname = []
        element = []
        charge = []
        moltype = []
        charmm_type = []
        atom_charge = []
        atom_vdw = []
        residue_flag = []
        original_index = []
        original_resid = []
        unique_names = []
        unique_resnames = []
        unique_resids = []
        unique_chains = []
        unique_segnames = []
        unique_occupancies = []
        unique_betas = []
        unique_elements = []
        unique_moltypes = []
        unique_charmm_type = []
        conect = {}
        natoms = self.natoms()

        for i in xrange(natoms):
            if(mask[i] == 1):
                try:
                    # if True:
                    atom.append(self._atom[i])
                    index.append(self._index[i])
                    original_index.append(self._original_index[i])
                    name.append(self._name[i])
                    loc.append(self._loc[i])
                    resname.append(self._resname[i])
                    chain.append(self._chain[i])
                    resid.append(self._resid[i])
                    original_resid.append(self._original_resid[i])
                    rescode.append(self._rescode[i])
                    occupancy.append(self._occupancy[i])
                    beta.append(self._beta[i])
                    segname.append(self._segname[i])
                    element.append(self._element[i])
                    charge.append(self._charge[i])
                    moltype.append(self._moltype[i])
                    residue_flag.append(self._residue_flag[i])
                    try:
                        charmm_type.append(self._charmm_type[i])
                    except:
                        pass
                    try:
                        atom_charge.append(self._atom_charge[i])
                    except:
                        pass
                    try:
                        atom_vdw.append(self._atom_vdw[i])
                    except:
                        pass

                    if(self._name[i] not in unique_names):
                        unique_names.append(self._name[i])
                    if(self._resname[i] not in unique_resnames):
                        unique_resnames.append(self._resname[i])
                    if(self._resid[i] not in unique_resids):
                        unique_resids.append(self._resid[i])
                    if(self._chain[i] not in unique_chains):
                        unique_chains.append(self._chain[i])
                    if(self._segname[i] not in unique_segnames):
                        unique_segnames.append(self._segname[i])
                    if(self._occupancy[i] not in unique_occupancies):
                        unique_occupancies.append(self._occupancy[i])
                    if(self._beta[i] not in unique_betas):
                        unique_betas.append(self._beta[i])
                    if(self._element[i] not in unique_elements):
                        unique_elements.append(self._element[i])
                    if(self._moltype[i] not in unique_moltypes):
                        unique_moltypes.append(self._moltype[i])

                except:
                    error.append(
                        'failed in copy_molecule when attempting to assign descriptors to atom ' + str(i))
                    print('\n\nerror = ', error)
                    sys.stdout.flush()
                    return error

        other.setAtom(atom)
        other.setIndex(numpy.array(index, numpy.int))

        original_index = numpy.array(original_index, numpy.int)
        other.setOriginal_index(original_index)

        other.setName(name)
        other.setLoc(loc)
        other.setResname(resname)
        other.setChain(chain)
        other.setResid(numpy.array(resid, numpy.int))
        other.setOriginal_resid(numpy.array(original_resid, numpy.int))
        other.setRescode(rescode)
        other.setOccupancy(occupancy)
        other.setBeta(beta)
        other.setSegname(segname)
        other.setElement(element)
        other.setCharge(charge)
        other.setMoltype(moltype)
        other.setCharmm_type(charmm_type)
        error, coor = self.get_coor_using_mask(frame, mask)
        other.setCoor(coor)
        other.setNatoms(len(index))
        other.setResidue_flag(residue_flag)

        other.setAtom_charge(numpy.array(atom_charge, numpy.float32))
        other.setAtom_vdw(numpy.array(atom_vdw, numpy.float32))

        other._number_of_names = len(unique_names)
        other._names = unique_names
        other._number_of_resnames = len(unique_resnames)
        other._resnames = unique_resnames
        other._number_of_resids = len(unique_resids)
        other._resids = unique_resids
        other._number_of_chains = len(unique_chains)
        other._chains = unique_chains
        other._number_of_segnames = len(unique_segnames)
        other._segnames = unique_segnames
        other._number_of_occupancies = len(unique_occupancies)
        other._occupancies = unique_occupancies
        other._number_of_betas = len(unique_betas)
        other._betas = unique_betas
        other._number_of_elements = len(unique_elements)
        other._elements = unique_elements
        other._number_of_moltypes = len(unique_moltypes)
        other._moltypes = unique_moltypes

        for oindex in original_index:

            if oindex in self.conect():

                conect[oindex] = self.conect()[oindex]

        other.setConect(conect)

        return error

    def duplicate_molecule(self, number_of_duplicates, **kwargs):
        '''
        This method copies all attributes from one molecule to a new
        set of a user-supplied number of duplicate molecules

        Parameters
        ----------
        number_of_duplicates 
            integer : number of copies to make

        kwargs 
            optional future arguments
                           
        Returns
        -------
        molecules
            list of system objects

        Examples
        --------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.coor()[0][0]
        array([-21.52499962, -67.56199646,  86.75900269])
        >>> molecule.name()[:10]
        ['N', 'HT1', 'HT2', 'HT3', 'CA', 'HA1', 'HA2', 'C', 'O', 'N']  
       
        >>> import sasmol.util as utilities
        >>> number_of_duplicates = 108
        >>> molecules = utilities.duplicate_molecule(molecule, number_of_duplicates)
        >>> molecules[-1].coor()[0][0]
        array([-21.52499962, -67.56199646,  86.75900269]) 
        >>> molecules[-1].name()[:10]
        ['N', 'HT1', 'HT2', 'HT3', 'CA', 'HA1', 'HA2', 'C', 'O', 'N']  
      
       
        Note
        ____
        Using deepcopy directly in subset.py leads to inheritance conflict.  
        Therefore subset calls a method held in utilities to make duplicates. 
        
        '''
    
        molecules = utilities.duplicate_molecule(molecule, number_of_duplicates)
         
        return molecules

    def get_indices_from_mask(self, mask):
        '''
        This method returns the internal indices for the supplied
        mask.  

        Parameters
        ----------
        mask 
            numpy integer array : mask array of length of the number of atoms
                                  with 1 or 0 for each atom depending on the selection
                                  used to create the mask

        kwargs 
            optional future arguments
                           
        Returns
        -------
        indices
            numpy integer array : indices of atoms determined by the input mask


        Examples
        --------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> basis_filter = "name[i] == 'CA'"
        >>> error, mask = molecule.get_subset_mask(basis_filter)
        >>> indices = molecule.get_indices_from_mask(mask) 
        >>> indices[:10]
        array([  4,  11,  21,  45,  55,  66,  82, 101, 112, 119])

        '''

        natoms = self.natoms()
        indices = numpy.nonzero(mask * numpy.arange(1, natoms + 1))[0]

        return indices

    def get_coor_using_mask(self, frame, mask):
        '''
        This method extracts coordinates from frame=frame of system object (self)
        using a supplied mask which has been created before this method is called.

        Coorindates are chosen for the elements that are equal to 1 in the supplied mask array.

        Parameters
        ----------
        frame 
            integer : trajectory frame number to use
        
        mask 
            numpy integer array : mask array of length of the number of atoms
                                  with 1 or 0 for each atom depending on the selection
                                  used to create the mask

        kwargs 
            optional future arguments
                           
        Returns
        -------
        error
            string : error statement

        coor
            coordinates corresponding to those determined by the input mask

        Examples
        --------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> basis_filter = "name[i] == 'CA'"
        >>> error, mask = molecule.get_subset_mask(basis_filter)
        >>> frame = 0
        >>> error, coor = molecule.get_coor_using_mask(frame, mask)
        >>> coor[0][0]
        array([-21.72500038, -66.91000366,  85.45700073], dtype=float32)

        '''
        error = []
        new_coor = []
        natoms = self.natoms()

        indicies = numpy.nonzero(mask * numpy.arange(1, natoms + 1))[0]

        this_frame_coor = self._coor[frame, :, :]

        coor = numpy.zeros((1, len(indicies), 3), numpy.float32)

        try:
            coor[0] = numpy.take(this_frame_coor[:, :], indicies, 0)
        except:
            error.append(
                'failed to extract coordinates from frame ' + str(frame))
            return error, coor

        return error, coor

    def set_coor_using_mask(self, other, frame, mask):
        '''
        This method replaces coordinates from frame=frame of system object (self)
        using a supplied mask which has been created before this method is called.

        Coordinates are chosen for the elements that are equal to 1 in the supplied mask array.

        Parameters
        ----------
        frame 
            integer : trajectory frame number to use
        
        mask 
            numpy integer array : mask array of length of the number of atoms
                                  with 1 or 0 for each atom depending on the selection
                                  used to create the mask

        kwargs 
            optional future arguments
                           
        Returns
        -------
        error
            string : error statement

            updated self._coor

        Examples
        --------

        >>> import sasmol.system as system
        >>> molecule_1 = system.Molecule('hiv1_gag.pdb')
        >>> molecule_2 = system.Molecule('other_hiv1_gag.pdb')
        >>> basis_filter = "name[i] == 'CA'"
        >>> error, mask = molecule_1.get_subset_mask(basis_filter)
        >>> frame = 0
        >>> error = molecule_1.set_coor_using_mask(molecule_2, frame, mask)

        Note
        ____
        molecule_2 must be smaller or equal to molecule_1 and that the coordinates 
        in molecule_2 are in the same order in molecule_1 

        '''
        error = []
        natoms_self = self.natoms()
        indicies_self = numpy.nonzero(
            mask * numpy.arange(1, natoms_self + 1))[0]

        three_indicies_self = []
        for i in xrange(len(indicies_self)):
            this_index = indicies_self[i] * 3
            three_indicies_self.append(
                [this_index, this_index + 1, this_index + 2])
        three_indicies_self = numpy.array(three_indicies_self).flatten()

        coor = self.coor()[:, :, :]

        try:
            numpy.put(coor[frame], three_indicies_self,
                      other._coor[frame, :, :])

            self.setCoor(coor)

        except:
            error.append(
                'failed to replace coordinates from frame ' + str(frame))
            return error

        return error

    def set_descriptor_using_mask(self, mask, descriptor, value):
        '''
        This method writes the "value" to the given descriptor to
        the elements that are equal to 1 in the supplied mask array.

        Parameters
        ----------
        mask 
            numpy integer array : mask array of length of the number of atoms
                                  with 1 or 0 for each atom depending on the selection
                                  used to create the mask

        descriptor 
            system property : a property defined in an instance of a system object

        value
            string : new value to apply to selection defined by mask

        kwargs 
            point = True : will translate to a fixed point
                           given by value variable                                                                             
        Returns
        -------
        None
            updated self._descriptor

        Examples
        --------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.beta()[:10]
        ['0.00', '0.00', '0.00', '0.00', '0.00', '0.00', '0.00', '0.00', '0.00', '0.00']

        >>> basis_filter = "name[i] == 'CA'"
        >>> error, mask = molecule.get_subset_mask(basis_filter)
        >>> descriptor = molecule.beta()
        >>> value = '1.00'
        >>> error = molecule.set_descriptor_using_mask(mask, descriptor, value)
        >>> descriptor[:10]
        ['0.00', '0.00', '0.00', '0.00', '1.00', '0.00', '0.00', '0.00', '0.00', '0.00']
    
        which can then be used to set the new values into the molecule 

        >>> molecule.setBeta(descriptor)
        >>> molecule.beta()[:10]
        ['0.00', '0.00', '0.00', '0.00', '1.00', '0.00', '0.00', '0.00', '0.00', '0.00']
       
        Note
        ____
        Coordinate arrays can not be manipulated by this method.

        TODO: If possible, get rid of loop

        '''
        error = []

        natoms = self.natoms()
        for i in xrange(natoms):
            if(mask[i] == 1):
                try:
                    descriptor[i] = value
                except:
                    error.append('failed to assign ' +
                                 str(value) + ' to atom ' + str(i))
                    return error
        return error

    def apply_biomt(self, frame, selection, U, M, **kwargs):
        """
        Apply biological unit transforms (BIOMT) to the coordinates of the
        chosen selection and frame.

        Information on BIOMT available at:
        http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies

        Parameters
        ----------
        frame   
            integer : frame number with coordinates to transform

        selection
            string  : selection string in standard SASMOL format  
                      specifying the coordinates to be transformed

        U
            numpy array : 3 x 3 rotation matrix

        M
            numpy array : 3 x 1 translation vector

        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        None
            updated self._coor 

        Examples
        -------
        
        Note
        ____

        TODO: add example 

        """

        # Get the coordinates for just the chosen frame and selection
        # as a masked numpy array
        coords = self.coor()
        err, mask = self.get_subset_mask(selection)
        err, sel_coords = self.get_coor_using_mask(frame, mask)

        # Transform coordinates
        new_coords = numpy.dot(U, sel_coords[frame].T).T
        new_coords = new_coords + M

        # Re-write edited coordinates into original array
        coords[frame] = new_coords

        return

    def copy_apply_biomt(self, other, frame, selection, U, M, **kwargs):
        """
        Copy selected atoms (with initial coordinates from the given frame)
        to new Molecule object (other) and apply transforms taken from biological
        unit (BIOMT) to the coordinates.

        Information on BIOMT available at:
        http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies

        Parameters
        ----------
        other
            system object : object to copy transformed information into
        
        frame   
            integer : frame number with coordinates to transform

        selection
            string  : selection string in standard SASMOL format  
                      specifying the coordinates to be transformed

        U
            numpy array : 3 x 3 rotation matrix

        M
            numpy array : 3 x 1 translation vector

        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        None
            updated self._coor 

        Examples
        -------

        Note
        ____

        TODO: add example 

        """

        # Copy selected atoms to new molecule (other)
        err, mask = self.get_subset_mask(selection)
        err = self.copy_molecule_using_mask(other, mask, frame)

        err = other.apply_biomt(0, selection, U, M)

        return


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
