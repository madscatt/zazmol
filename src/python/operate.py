from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

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
        and alignment operatioins on instances of system objects.

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
        mass_check determines if mass is defined for the ojbect so that
        center of mass (COM) can be calculated
        

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
            print(type(self._com))

            self._coor[frame, :, 0] = self._coor[frame, :, 0] - self._com[0]
            self._coor[frame, :, 1] = self._coor[frame, :, 1] - self._com[1]
            self._coor[frame, :, 2] = self._coor[frame, :, 2] - self._com[2]
         
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
        >>> displacement = [3.0, 4.0, 5.0]
        >>> molecule.translate(frame,displacement)
        >>> molecule.calculate_center_of_mass(frame)
        array([ -3.79114736, -19.71577133,  13.06558513])

        Note 
        ----------
        mass_check is called to validate self._total_mass()    
        
        '''

        self.mass_check()
        self._com = self.calculate_center_of_mass(frame)

        self._coor[frame, :, 0] = self._coor[frame, :, 0] - self._com[0]
        self._coor[frame, :, 1] = self._coor[frame, :, 1] - self._com[1]
        self._coor[frame, :, 2] = self._coor[frame, :, 2] - self._com[2]

        self._com = self.calculate_center_of_mass(frame)

        return

    def align(self, frame, coor_sub_2, com_sub_2, coor_sub_1, com_sub_1, **kwargs):
        '''
            Alignment of one object on top of another
            "self" is aligned onto "other" using the basis
            of molecule 2 to align onto the basis of molecule 1
            and the transformation is then done to all the atoms of
            molecule 2

        '''
        self.mass_check()
        self.calculate_center_of_mass(frame)

        u = linear_algebra.find_u(coor_sub_1, coor_sub_2)

        tao = numpy.transpose(self.coor()[frame] - com_sub_2)

        error, nat2 = linear_algebra.matrix_multiply(u, tao)

        ncoor = numpy.transpose(nat2) + com_sub_1

        self._coor[frame, :] = ncoor

        return

    def rotate(self, frame, axis, theta, **kwargs):
        '''
            Simple rotation about the x, y, or z axis.

            Note that calcuations are in radians

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

    def general_axis_rotate(self, frame, theta, ux, uy, uz, **kwargs):
        '''
            The general rotation of a molecule along an arbitrarily
            given unit axis (ux,uy,uz) by an angle theta.

        Note that calcuations are in radians
        '''

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

    def euler_rotate(self, frame, phi, theta, psi, **kwargs):
        '''
        Rotate the molecule by a euler angle set (phi,theta,psi)

        Note that calcuations are in radians

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
