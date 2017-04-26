from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
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
#	CALCULATE
#
#	12/10/2009	--	initial coding			        :	jc
#	12/20/2015	--	refactored for release          :	jc
#	07/23/2016	--	refactored for Python 3         :	jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **

'''
    Calculate contains the classes and methods to calculate various
    atomic and molecular properties from instances of system objects

'''

import sys
import numpy
import sasmol.operate as operate

class Calculate(object):

    """ Base class containing methods to calculate properties of system object.

        Examples
        ========

        First example shows how to use class methods from system object:

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_mass()
        47896.61864599498

        Second example shows how to use class methods directly:

        >>> import sasmol.system as system
        >>> import sasmol.calculate as calculate
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> calculate.Calculate.calculate_mass(molecule)
        47896.61864599498

        Note
        ----

        `self` parameter is not shown in the ``Parameters`` section in the documentation

        TODO:  Need to write a generic driver to loop over single or multiple frames

    """

    def calculate_mass(self, **kwargs):
        '''

        Note
        ----

        atomic weights are contained in the ``properties.py`` file

        http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some

        standard atomic weight is based on the natural istopic composition

        NOTE: deuterium is 'D' 2H1 and '1H' is 1H1, all other elements
        have their natural abundance weight. These elements are located
        at the end of the dictionary.

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        float
            mass of object in Daltons

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_mass()
        47896.61864599498

        '''

        standard_atomic_weight = self.amu()
        self._total_mass = 0.0
        self._mass = numpy.zeros(len(self._element), numpy.float)

        count = 0

        for element in self._element:
            if element in standard_atomic_weight:
                self._total_mass = self._total_mass + \
                    standard_atomic_weight[element]
                self._mass[count] = standard_atomic_weight[element]
                count += 1
            else:
                message = 'element ' + element + ' not found'
#                log.error('ERROR: ' + message)

#               need to return an error that element was not found

        return self._total_mass

    def calculate_center_of_mass(self, frame, **kwargs):
        '''
        This method calculates the center of mass of the object.

        Parameters
        ----------
        frame
            integer : trajectory frame number to use

        kwargs
            optional future arguments

        Returns
        -------
        numpy array
            coordinates of center of mass

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_center_of_mass(0)
        array([ -6.79114736, -23.71577133,   8.06558513])

        '''

        if(self._total_mass <= 0.0):
            self.calculate_mass()

        x = self._coor[frame, :, 0]
        y = self._coor[frame, :, 1]
        z = self._coor[frame, :, 2]

        comx = numpy.sum(self._mass * x) / self._total_mass
        comy = numpy.sum(self._mass * y) / self._total_mass
        comz = numpy.sum(self._mass * z) / self._total_mass

        self._com = numpy.array([comx, comy, comz], numpy.float)

        return self._com

    def calculate_radius_of_gyration(self, frame, **kwargs):
        '''
        This method calculates the radius of gyration of the object

        Parameters
        ----------
        frame
            integer : trajectory frame number to use

        kwargs
            optional future arguments

        Returns
        -------
        float
            radius of gyration of object

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_radius_of_gyration(0)
        64.043168998442368

        '''

        self._com = self.calculate_center_of_mass(frame)

        if(self._natoms > 0):
            rg2 = ((self._coor[frame, :, :] - self._com)
                    * (self._coor[frame, :, :] - self._com))
            self._rg = numpy.sqrt(numpy.sum(numpy.sum(rg2)) / self._natoms)

        return self._rg

    def calculate_root_mean_square_deviation(self, other, **kwargs):
        '''
        This method calculates the radius root mean square deviation (rmsd)
        of one set of coordinates compared to another

        self contains the coordinates of set 1
        other contains the coordinates of set 2

        the number of coordinates in each set must be equal

        To use this over multiple frames you must call this function
        repeatedly.

        Parameters
        ----------
        other
            system object with coordinates with equal number of frames

        kwargs
            optional future arguments

        Returns
        -------
        float
            root mean square deviation between objects

        Examples
        -------


        '''

        # OPEN	Add frame support here?
        try:
            dxyz = ((self._coor - other._coor) * (self._coor - other._coor))
            self._rmsd = numpy.sqrt((numpy.sum(dxyz)) / self._natoms)
        except:
            if(self._natoms != other._natoms):
                print('number of atoms in (1) != (2)')
                print('rmsd not calculated: None returned')
                print('number of atoms in self is < 1')
                print('number of atoms in other is < 1')
                self._rmsd = None

        return self._rmsd

    def calculate_principal_moments_of_inertia(self, frame, **kwargs):
        '''
        This method calculates the principal moments of inertia
        of the object. It uses the center method from operate.Move
        to center the molecule.
        The present method is designated for the non-linear system with
        non-singular moment of inertia matrix only. For the linear systems, it
        will return eigenvectors and I as None.  Testing for non-None return
        values should be done in the calling method.

        Parameters
        ----------
        frame
            integer : trajectory frame number to use

        kwargs
            optional future arguments

        Returns
        -------
        tuple of numpy arrays
            principal moments of inertia of object :
            eigenvalues, eigenvectors, and I

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_principal_moments_of_inetia(0)
        (array([  1.30834716e+07,   1.91993314e+08,   1.85015201e+08]),
        array([[-0.08711655, -0.97104917,  0.22242802],
               [-0.512547  ,  0.23514759,  0.82583363],
               [ 0.85422847,  0.04206103,  0.51819358]]),
        array([[  1.90290278e+08,  -9.27036144e+06,   1.25097100e+07],
               [ -9.27036144e+06,   1.40233826e+08,   7.53462715e+07],
               [  1.25097100e+07,   7.53462715e+07,   5.95678834e+07]]))

        '''

        com = self.calculate_center_of_mass(frame)

        operate.Move.center(self, frame)

        n_atoms = self._natoms
        m = self._mass.reshape(n_atoms, -1)
        m_coor = m * self._coor[frame]
        m_coor2 = numpy.dot(self._coor[frame].T, m_coor)
        numpy.fill_diagonal(m_coor2, m_coor2.diagonal() - m_coor2.trace())
        I = -m_coor2

        if numpy.linalg.matrix_rank(I) < 3:
            print("You are requesting the pmi calculation for a singular system.")
            print("The eigen-decomposition of this system is not defined")

            uk = None
            ak = None
            I = None
        else:
            uk, ak = numpy.linalg.eig(I)
            order = uk.argsort()
            uk = uk[order]
            ak = ak[:, order]

        operate.Move.translate(self, frame, com, point=True)

        return uk, ak, I

    def calculate_minimum_and_maximum(self, **kwargs):
        '''
        This method calculates the minimum and maximum values
        of of the object in (x,y,z)

        The default usage is to evaluate all frames

        A numpy array of minimum and maximum values for each
        dimension are returned

        Parameters
        ----------
        kwargs
           frames = [] : integer list of frames to process

        Returns
        -------
        numpy array
            nested list of minimum and maximum values
            [ [ min_x, min_y, min_z ], [max_x, max_y, max_z] ]

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_minimum_and_maximum()
        [array([-31.29899979, -93.23899841, -85.81900024]), array([ 19.64699936,  30.37800026,  99.52999878])]

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.read_dcd('hiv1_gag_200_frames')
        >>> molecule.calculate_minimum_and_maximum()
        [array([ -94.47146606, -121.88082886,  -99.94940948]), array([  52.85133362,   65.53725433,  100.76850891])]
        >>> molecule.calculate_minimum_and_maximum(frames=[0,1,2,3])
        [array([-30.9330883 , -92.68256378, -84.51082611]), array([ 20.98281288,  38.45230484,  99.91564178])]


        '''

        try:
            frames = kwargs['frames']
        except:
            frames = [x for x in xrange(self.number_of_frames())]

        first_flag = True

        for frame in frames:

            this_min_x=numpy.min(self._coor[frame,:,0])
            this_max_x=numpy.max(self._coor[frame,:,0])
            this_min_y=numpy.min(self._coor[frame,:,1])
            this_max_y=numpy.max(self._coor[frame,:,1])
            this_min_z=numpy.min(self._coor[frame,:,2])
            this_max_z=numpy.max(self._coor[frame,:,2])

            if(first_flag or (this_min_x < min_x)):
                min_x = this_min_x
            if(first_flag or (this_min_y < min_y)):
                min_y = this_min_y
            if(first_flag or (this_min_z < min_z)):
                min_z = this_min_z

            if(first_flag or (this_max_x > max_x)):
                max_x = this_max_x
            if(first_flag or (this_max_y > max_y)):
                max_y = this_max_y
            if(first_flag or (this_max_z > max_z)):
                max_z = this_max_z

            first_flag = False

        self._minimum = numpy.array([min_x, min_y, min_z])
        self._maximum = numpy.array([max_x, max_y, max_z])

        return [self._minimum, self._maximum]

    def calculate_residue_charge(self, **kwargs):
        '''

        Method to sum the atomic charges and assign the net charge of the
        resiude to a new variable that is attached to each atom.

        Note
        ----------
        Requires that the atom_charge() attribute of object is complete

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        float
            charge per residue assigned to object

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_residue_charge()
        >>> single_molecule.calculate_residue_charge()
        >>> residue_charge = single_molecule.residue_charge()
        >>> print('res-charge = ',residue_charge[0])

        '''

        resid = self.resid()
        natoms = self.natoms()

        ###
        # needs a gentle failure if atom_charge() does not exist
        ###

        atom_charge = self.atom_charge()

        charge_sum = 0.0
        charge_residue_sum = []
        last_resid = resid[0]

        for i in xrange(natoms):
            this_resid = resid[i]
            this_charge = atom_charge[i]

            if(this_resid != last_resid or i == natoms - 1):
                charge_residue_sum.append([last_resid, charge_sum])
                charge_sum = this_charge
                last_resid = this_resid
            else:
                charge_sum += this_charge

        last_resid = resid[0]
        charge_residue = []

        for i in xrange(natoms):
            this_resid = resid[i]
            for j in xrange(len(charge_residue_sum)):
                if(this_resid == charge_residue_sum[j][0]):
                    charge_residue.append(charge_residue_sum[j][1])
                    continue

        self.setResidue_charge(numpy.array(charge_residue, numpy.float32))

        return

    def calculate_molecular_formula(self, **kwargs):
        '''

        Method to determine the number of each element in the molecule

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        dictionary
            {element : integer number, ... }

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.calculate_molecular_formula()
        {'H': 3378, 'C': 2080, 'S': 24, 'O': 632, 'N': 616}

        '''

        my_formula = {}

        for element in self._element:

            if element in my_formula:
                my_formula[element] += 1
            else:
                my_formula[element] = 1

        self.setFormula(my_formula)

        return my_formula
