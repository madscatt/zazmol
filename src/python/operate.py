from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

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
import sys
import numpy
import sasmol as sasmol
import sasmol.linear_algebra as linear_algebra

#	OPERATE
#
#	12/13/2009	--	initial coding			        :	jc
#   07/23/2016  --  refactored for release          :   jc
#   07/23/2016  --  refactored for Python 3         :   jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Operate contains the classes and methods to perform the basic
	translation, rotation, and alignment operations on instances of objects 
	described in the sasmol module.

	The methods in this class move entire objects, not pieces
	of single objects (i.e. no intra-object movements)
'''


class Move():

    """ Base class containing methods to perform basic translation, rotation,
        and alignment operations on instances of system objects.

        The methods in this class move entire objects, not pieces
        of single objects (i.e. no intra-object movements)

        Examples
        ========

        First example shows how to use class methods from system object:

        >>> import sasmol.system as system ; import math
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0 ; axis = 'x' ; theta = 45.0 * math.pi / 180.0
        >>> molecule.rotate(frame, axis, theta)

        Second example shows how to use class methods directly:

        >>> import sasmol.system as system ; import math
        >>> import sasmol.operate as operate
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0 ; axis = 'x' ; theta = 45.0 * math.pi / 180.0
        >>> operate.Move.rotate(molecule, frame, axis, theta) 

        Note
        ----

        `self` parameter is not shown in the ``Parameters`` section in the documentation


    """

    '''


        moveto moves the COM to a point in space (x,y,z)

    '''

    def mass_check(self, **kwargs):
        ''' 
        Note
        ----
        mass_check determines if mass is defined for the object so that
        center of mass can be calculated
        

        Parameters
        ----------
        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        None
            updated self._total_mass

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.mass_check()
        
        ''' 
          
        if(self._total_mass <= 0.0):
            self.calculate_mass()
        return

    def translate(self, frame, value, **kwargs):
        ''' 
        
        translate moves the object 
        
        Parameters
        ----------
        frame 
            integer : trajectory frame number to use

        value 
            list of x, y, z float values

        kwargs 
            point = True : will translate to a fixed point
                           given by value variable                                                                             
        Returns
        -------
        None
            updated self._coor and self._center_of_mass

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0
        >>> molecule.calculate_center_of_mass(frame)
        array([ -6.79114736, -23.71577133,   8.06558513])
        >>> displacement = [3.0, 4.0, 5.0]
        >>> molecule.translate(frame, displacement)
        >>> molecule.calculate_center_of_mass(frame)
        array([ -3.79114736, -19.71577133,  13.06558513])
      
         
        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0
        >>> final_position = [3.0, 4.0, 5.0]
        >>> molecule.translate(frame, final_position, point=True)
        >>> molecule.calculate_center_of_mass(frame)  
        array([ 3.,  4.,  5.])
         
        Note 
        ----------
        mass_check is called to validate self._total_mass()    


        '''
        try: 
            point_flag = kwargs['point']
        except:
            point_flag = False
        
        self.mass_check()

        if point_flag:
            self._com = self.calculate_center_of_mass(frame)

            self._coor[frame, :, 0] -= self._com[0]
            self._coor[frame, :, 1] -= self._com[1]
            self._coor[frame, :, 2] -= self._com[2]
         
        self._coor[frame, :, 0] += value[0]
        self._coor[frame, :, 1] += value[1]
        self._coor[frame, :, 2] += value[2]

        self._com = self.calculate_center_of_mass(frame)

        return

    def center(self, frame, **kwargs):
        '''
        Method moves the center of mass of object to [0.0, 0.0, 0.0]
       
        Parameters
        ----------
        frame 
            integer : trajectory frame number to use

        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        None
            updated self._coor and self._center_of_mass

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0
        >>> molecule.calculate_center_of_mass(frame)
        array([ -6.79114736, -23.71577133,   8.06558513])
        >>> molecule.center(frame)
        >>> molecule.calculate_center_of_mass(frame)
        array([  7.11544707e-13,   2.48159571e-12,  -8.45832820e-13])

        Note 
        ----------
        mass_check is called to validate self._total_mass()    
        
        Can achieve same result using self.translate(frame, [0,0,0], point=True)

        
        '''

        self.mass_check()
        self._com = self.calculate_center_of_mass(frame)

        self._coor[frame, :, 0] = self._coor[frame, :, 0] - self._com[0]
        self._coor[frame, :, 1] = self._coor[frame, :, 1] - self._com[1]
        self._coor[frame, :, 2] = self._coor[frame, :, 2] - self._com[2]

        self._com = self.calculate_center_of_mass(frame)

        return
    
    def align(self, other, self_basis, other_basis, **kwargs):
        '''
            Alignment of one object on top of another
            "self" is aligned onto "other" using the basis
            of molecule 2 to align onto the basis of molecule 1
            and the transformation is then done to all the atoms of
            molecule 2

            self = molecule_2

            other = molecule_1

            self aligned to other

            molecule_2 aligned to molecule_1


        Parameters
        ----------
        frame 
            integer : trajectory frame number to use
        
        other
            system object : molecule 1

        self_basis
            string : unique description of atoms used for alignment

        other_basis
            string : unique description of atoms used for alignment
        
        kwargs 
            optional future arguments

        Returns
        -------
        None
            updated self._coor

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule_1 = system.Molecule('hiv1_gag.pdb')
        >>> molecule_2 = system.Molecule('moved_and_rotated_hiv1_gag.pdb')
        >>> frame = 0
        >>> basis_1 = 'name[i] == "CA"'
        >>> basis_2 = 'name[i] == "CA"'
        >>> molecule_2.align(molecule_1, basis_1, basis_2)
        >>> com_sub_2 = molecule_2.calculate_center_of_mass(frame)
        
        Note
        ----
        mass_check determines if mass is defined for the object so that
        center of mass can be calculated
         
        
        '''

        frame = 0


        ### other = molecule_1 (reference)
        
        error, other_mask = other.get_subset_mask(other_basis)
        
        subset_other = sasmol.system.Molecule()
        error = other.copy_molecule_using_mask(subset_other, other_mask, frame) 
        
        com_subset_other = subset_other.calculate_center_of_mass(frame)
        subset_other.center(frame)
        coor_subset_other = subset_other.coor()[frame]
       

        ### self = molecule_2 (to be aligned to other / molecule_1)

        error, self_mask = self.get_subset_mask(self_basis)
        
        subset_self = sasmol.system.Molecule()
        error = self.copy_molecule_using_mask(subset_self, self_mask, frame) 
        
        com_subset_self = subset_self.calculate_center_of_mass(frame)
        subset_self.center(frame)
        coor_subset_self = subset_self.coor()[frame]
       
        
        u = linear_algebra.find_u(coor_subset_self, coor_subset_other)

        tao = numpy.transpose(self.coor()[frame] - com_subset_other)

        error, nat2 = linear_algebra.matrix_multiply(u, tao)

        ncoor = numpy.transpose(nat2) + com_subset_other

        self._coor[frame, :] = ncoor

        return

    def rotate(self, frame, axis, theta, **kwargs):
        '''
            Simple rotation about the x, y, or z axis.


        Parameters
        ----------
        frame 
            integer : trajectory frame number to use

        axis 
            string : 'x', 'y', or 'z'

        theta 
            float : angle in radians

        kwargs 
            optional future arguments        
       
        Returns
        -------
        None
            updated self._coor 

        Examples
        -------

        >>> import sasmol.system as system ; import math
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0 ; axis = 'x' ; theta = 45.0 * math.pi / 180.0
        >>> molecule.rotate(frame, axis, theta)
 
        Note 
        ----------
        Calculations are carried out using radians
        
        '''

        cs = numpy.cos(theta)
        si = numpy.sin(theta)
        if(axis == 'x'):
            mat = numpy.array(
                    [[1.0, 0.0, 0.0], [0.0, cs, -si], [0.0, si, cs]])
        elif(axis == 'y'):
            mat = numpy.array(
                    [[cs, 0.0, si], [0.0, 1.0, 0.0], [-si, 0.0, cs]])
        elif(axis == 'z'):
            mat = numpy.array(
                    [[cs, -si, 0.0], [si, cs, 0.0], [0.0, 0.0, 1.0]])

        coordt = self._coor[frame, :].T
        error, matrix_product = linear_algebra.matrix_multiply(mat, coordt)
        ncoord = matrix_product.T
        self._coor[frame, :] = ncoord

        return

    def rotate_general_axis(self, frame, theta, unit_axis, **kwargs):
        '''
            The general rotation of a molecule along an arbitrarily
            given unit axis (ux,uy,uz) by an angle theta.

        Parameters
        ----------
        frame 
            integer : trajectory frame number to use

        theta 
            float : angle in radians
        
        unit_axis
            float : [ux, uy, uz] components of unit axis

        kwargs 
            optional future arguments        
       
        Returns
        -------
        None
            updated self._coor 

        Examples
        -------

        >>> import sasmol.system as system ; import math
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0 ; theta = 45.0 * math.pi / 180.0
        >>> unit_axis = [ 0.2, 1.3, -3.5 ]
        >>> molecule.rotate_general_axis(frame, theta, unit_axis)
 
        Note 
        ----------
        Calculations are carried out using radians
        
        '''

        ux = unit_axis[0]
        uy = unit_axis[1]
        uz = unit_axis[2]

        c11 = numpy.cos(theta) + pow(ux, 2) * (1 - numpy.cos(theta))
        c12 = ux * uy * (1 - numpy.cos(theta)) - uz * numpy.sin(theta)
        c13 = ux * uz * (1 - numpy.cos(theta)) + uy * numpy.sin(theta)
        c21 = uy * ux * (1 - numpy.cos(theta)) + uz * numpy.sin(theta)
        c22 = numpy.cos(theta) + pow(uy, 2) * (1 - numpy.cos(theta))
        c23 = uy * uz * (1 - numpy.cos(theta)) - ux * numpy.sin(theta)
        c31 = uz * ux * (1 - numpy.cos(theta)) - uy * numpy.sin(theta)
        c32 = uz * uy * (1 - numpy.cos(theta)) + ux * numpy.sin(theta)
        c33 = numpy.cos(theta) + pow(uz, 2) * (1 - numpy.cos(theta))

        C = numpy.matrix(
            [[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]])

        coor = numpy.array(self.coor()[frame] * C)

        self.coor()[frame, :] = coor

        return

    def rotate_euler(self, frame, phi, theta, psi, **kwargs):
        '''
        Rotate the molecule by a euler angle set (phi,theta,psi)

        Parameters
        ----------
        frame 
            integer : trajectory frame number to use

        phi 
            float : phi angle

        theta 
            float : theta angle

        psi 
            float : psi angle

        kwargs 
            optional future arguments        
       
        Returns
        -------
        None
            updated self._coor 

        Examples
        -------

        >>> import sasmol.system as system ; import math
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> frame = 0 
        >>> phi = 45.0 * math.pi / 180.0 
        >>> theta = 45.0 * math.pi / 180.0
        >>> psi = 32.0 * math.pi / 180.0
        >>> molecule.euler_rotate(frame, phi, theta, psi)
 
        Note 
        ----------
        Calculations are carried out using radians

        '''

        c11 = numpy.cos(theta) * numpy.cos(psi)
        c12 = numpy.cos(phi) * numpy.sin(psi) + \
                numpy.sin(phi) * numpy.sin(theta) * numpy.cos(psi)
        c13 = numpy.sin(phi) * numpy.sin(psi) - \
                numpy.cos(phi) * numpy.sin(theta) * numpy.cos(psi)
        c21 = -numpy.cos(theta) * numpy.sin(psi)
        c22 = numpy.cos(phi) * numpy.cos(psi) - \
                numpy.sin(phi) * numpy.sin(theta) * numpy.sin(psi)
        c23 = numpy.sin(phi) * numpy.cos(psi) + \
                numpy.cos(phi) * numpy.sin(theta) * numpy.sin(psi)
        c31 = numpy.sin(theta)
        c32 = -numpy.sin(phi) * numpy.cos(theta)
        c33 = numpy.cos(phi) * numpy.cos(theta)

        C = numpy.matrix(
                [[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]])

        coor = numpy.array(self.coor()[frame] * C)

        self.coor()[frame, :] = coor

        return
