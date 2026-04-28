# from __future__ import unicode_literals

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
import math
import sasmol as sasmol
import sasmol.config as config
import sasmol.linear_algebra as linear_algebra

# OPERATE
#
# 12/13/2009	--	initial coding			        :	jc
#   07/23/2016  --  refactored for release          :   jc
#   07/23/2016  --  refactored for Python 3         :   jc
#
# 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
# *      **
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

        if (self._total_mass <= 0.0):
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
        except Exception:
            point_flag = False

        self.mass_check()

        if point_flag:
            self._center_of_mass = self.calculate_center_of_mass(frame)

            self._coor[frame, :, 0] -= self._center_of_mass[0]
            self._coor[frame, :, 1] -= self._center_of_mass[1]
            self._coor[frame, :, 2] -= self._center_of_mass[2]

        self._coor[frame, :, 0] += value[0]
        self._coor[frame, :, 1] += value[1]
        self._coor[frame, :, 2] += value[2]

        self._center_of_mass = self.calculate_center_of_mass(frame)

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
        self._center_of_mass = self.calculate_center_of_mass(frame)

        self._coor[frame, :, 0] = self._coor[frame, :, 0] - \
            self._center_of_mass[0]
        self._coor[frame, :, 1] = self._coor[frame, :, 1] - \
            self._center_of_mass[1]
        self._coor[frame, :, 2] = self._coor[frame, :, 2] - \
            self._center_of_mass[2]

        self._center_of_mass = self.calculate_center_of_mass(frame)

        return

    def align(self, other, self_basis=None, other_basis=None, mode='production', align_variables=None, **kwargs):
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
                string : unique description of atoms used for alignment (only needed in initialization mode)

            other_basis
                string : unique description of atoms used for alignment (only needed in initialization mode)

            mode
                string : 'initialization' or 'production'

            align_variables
                dict : data from initialization mode to be used in production mode

            kwargs 
                optional future arguments

            Returns
            -------
            None or dict
                updated self._coor or initialization data

            Examples
            -------

            >>> import sasmol.system as system
            >>> molecule_1 = system.Molecule('hiv1_gag.pdb')
            >>> molecule_2 = system.Molecule('moved_and_rotated_hiv1_gag.pdb')
            >>> frame = 0
            >>> basis_1 = 'name[i] == "CA"'
            >>> basis_2 = 'name[i] == "CA"'
            >>> align_variables = molecule_2.align(molecule_1, basis_1, basis_2, mode='initialization')
            >>> molecule_2.align(molecule_1, align_variables=align_variables)
            >>> com_sub_2 = molecule_2.calculate_center_of_mass(frame)

            Note
            ----
            mass_check determines if mass is defined for the object so that
            center of mass can be calculated
            '''
        frame = kwargs.get('frame', 0)

        if mode == 'initialization':
            if self_basis is None or other_basis is None:
                raise ValueError(
                    "self_basis and other_basis must be provided in initialization mode")

            # other = molecule_1 (reference)
            error, other_mask = other.get_subset_mask(other_basis)
            subset_other = sasmol.system.Molecule()
            error = other.copy_molecule_using_mask(
                subset_other, other_mask, frame)
            com_subset_other = subset_other.calculate_center_of_mass(frame)
            subset_other.center(frame)
            coor_subset_other = subset_other.coor()[frame]

            # self = molecule_2 (to be aligned to other / molecule_1)
            error, self_mask = self.get_subset_mask(self_basis)
            subset_self = sasmol.system.Molecule()
            error = self.copy_molecule_using_mask(
                subset_self, self_mask, frame)
            com_subset_self = subset_self.calculate_center_of_mass(frame)
            subset_self.center(frame)
            coor_subset_self = subset_self.coor()[frame]

            # Return initialization data as align_variables
            align_variables = {
                'self_mask': self_mask,
                'subset_self': subset_self,
                'coor_subset_other': coor_subset_other,
                'com_subset_other': com_subset_other,
            }

            return align_variables

        else:  # Default to production mode
            if align_variables is None:
                raise ValueError(
                    "align_variables must be provided in production mode")

            # assign initialization data to variables

            #frame = 0
            #added for sassie testing
            frame = kwargs.get('frame', 0)

            com_subset_other = align_variables['com_subset_other']
            coor_subset_other = align_variables['coor_subset_other']

            self_mask = align_variables['self_mask']
            subset_self = align_variables['subset_self']

            # update self object to current state
            error = self.copy_molecule_using_mask(
                subset_self, self_mask, frame)
            com_subset_self = subset_self.calculate_center_of_mass(frame)
            subset_self.center(frame)
            coor_subset_self = subset_self.coor()[frame]

            # u = linear_algebra.find_u(coor_subset_self, coor_subset_other)
            u = linear_algebra.find_u(coor_subset_other, coor_subset_self)
            tao = numpy.transpose(self.coor()[frame] - com_subset_self)
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
        if (axis == 'x'):
            mat = numpy.array(
                [[1.0, 0.0, 0.0], [0.0, cs, -si], [0.0, si, cs]])
        elif (axis == 'y'):
            mat = numpy.array(
                [[cs, 0.0, si], [0.0, 1.0, 0.0], [-si, 0.0, cs]])
        elif (axis == 'z'):
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

        C = numpy.array(
            [[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]],
            config.CALC_DTYPE)

        coor = numpy.dot(self.coor()[frame], C)

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

        C = numpy.array(
            [[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]],
            config.CALC_DTYPE)

        coor = numpy.dot(self.coor()[frame], C)

        self.coor()[frame, :] = coor

        return

    def align_pmi_on_cardinal_axes(self, frame):
        '''
        aligns all principal moment eigenvectors to cardinal x, y, & z axes
        '''

        self.align_pmi_on_axis(frame, 2, 'z')
        self.align_pmi_on_axis(frame, 1, 'y')

        return

    def align_pmi_on_axis(self, frame, pmi_eigenvector, alignment_vector_axis):
        '''
            aligns one principal moment eigenvector to the given axis
        '''
        self.center(frame)
        uk, ak, I = self.calculate_principal_moments_of_inertia(frame)

        ''' right handed coordinate frame '''
        ak = ak[:, pmi_eigenvector]

        axes = {
            'x': numpy.array([1.0, 0.0, 0.0]),
            'y': numpy.array([0.0, 1.0, 0.0]),
            'z': numpy.array([0.0, 0.0, 1.0]),
        }
        axis = axes[alignment_vector_axis]

        ak = ak / numpy.linalg.norm(ak)
        rotvec = numpy.cross(axis, ak)
        sine = numpy.linalg.norm(rotvec)
        cosine = numpy.clip(numpy.dot(ak, axis), -1.0, 1.0)

        if sine > 1.0e-12:
            rotvec = rotvec / sine
        elif cosine > 0.0:
            return
        else:
            basis = numpy.eye(3)[numpy.argmin(numpy.abs(ak))]
            rotvec = numpy.cross(ak, basis)
            rotvec = rotvec / numpy.linalg.norm(rotvec)

        # theta = math.atan(sine / cosine)
        theta = math.atan2(sine, cosine)

        unit_axis = [rotvec[0], rotvec[1], rotvec[2]]
        self.rotate_general_axis(frame, theta, unit_axis)

        return


def set_average_vdw(mol):
    '''
    # NOTE: atom_vdw must be shape (natoms,2); only column 2 (radius) is used by pairs.f
    '''

    element_vdw = []

    for e in mol.element():
        item = legacy_van_der_waals.get(e, None)

        if isinstance(item, list):
            element_vdw.append(item)
        else:
            element_vdw.append([None, None])

    mol.setAtom_vdw(element_vdw)

    return

# ----------------------------------------------------------------------
# LEGACY van_der_waals DATA (DO NOT DELETE WITHOUT REVIEW)
#
# This is a mixed legacy dataset historically used by set_average_vdw.
#
# IMPORTANT:
# - Some entries are LISTS: [A, B]
#     These likely represent Lennard-Jones–type parameters
#     (e.g., energy term and distance parameter, or A/B coefficients).
#
# - Some entries are SCALARS:
#     These are NOT vdW radii and appear to be atomic weights or placeholders.
#     They were effectively ignored by the original implementation.
#
# - Original behavior:
#     Only list entries (len(item) > 1) were used to populate atom_vdw.
#     Scalar entries were skipped or resulted in [None, None].
#
# - This dataset is NOT equivalent to:
#     properties.Atomic().van_der_waals_radii()
#
# - This table should be preserved for:
#     backward compatibility, historical reproducibility, and
#     potential force-field parameter interpretation.
#
# TODO (future, not now):
# - Determine exact physical meaning of [A, B]
# - Separate into:
#     * true vdW radii
#     * force-field parameters (ε/σ or A/B)
# ----------------------------------------------------------------------

legacy_van_der_waals = {
    'H': [-0.0368384615385, 0.928619230769],
    'He': 4.0026022,
    'Li': 6.9412, 'Be': 9.0121823, 'B': 10.8117,
    'C': [-0.0763111111111, 2.00249333333],
    'N': [-0.2, 1.85],
    'O': [-0.13805625, 1.7392625],
    'F': [-0.105, 1.7],
    'Ne': [-0.086000, 1.5300],
    'Na': [-0.0469, 1.36375],
    'Mg': [-0.0150, 1.18500],
    'Al': 26.9815382, 'Si': 28.08553,
    'P': [-0.585, 2.15],
    'S': [-0.450000, 2.000000],
    'Cl': [-0.150, 2.27],
    'Ar': 39.9481,
    'K': [-0.0870, 1.76375],
    'Ca': [-0.120, 1.367],
    'Sc': 44.9559108, 'Ti': 47.8671,
    'V': 50.94151, 'Cr': 51.99616, 'Mn': 54.9380499,
    'Fe': [0.000000, 0.650000],
    'Co': 58.9332009, 'Ni': 58.69342, 'Cu': 63.5463,
    'Zn': [-0.250000, 1.09000],
    'Ga': 69.7231, 'Ge': 72.641, 'As': 74.921602, 'Se': 78.963,
    'Br': 79.9041, 'Kr': 83.7982, 'Rb': 85.46783, 'Sr': 87.621,
    'Y': 88.905852, 'Zr': 91.2242, 'Nb': 92.906382, 'Mo': 95.942,
    'Tc': 98.0, 'Ru': 101.072, 'Rh': 102.905502, 'Pd': 106.421,
    'Ag': 107.86822, 'Cd': 112.4118, 'In': 114.8183, 'Sn': 118.7107,
    'Sb': 121.7601, 'Te': 127.603, 'I': 126.904473, 'Xe': 131.2936,
    'Cs': [-0.1900, 2.100],
    'Ba': 137.3277, 'La': 138.90552, 'Ce': 140.1161,
    'Pr': 140.907652, 'Nd': 144.243, 'Pm': 145.0, 'Sm': 150.363,
    'Eu': 151.9641, 'Gd': 157.253, 'Tb': 158.925342, 'Dy': 162.5001,
    'Ho': 164.930322, 'Er': 167.2593, 'Tm': 168.93421, 'Yb': 173.043,
    'Lu': 174.9671, 'Hf': 174.9671, 'Ta': 180.94791, 'W': 183.841,
    'Re': 186.2071, 'Os': 190.233, 'Ir': 192.2173, 'Pt': 195.0782,
    'Au': 196.966552, 'Hg': 200.592, 'Tl': 204.38332, 'Pb': 207.21,
    'Bi': 208.980382, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0, 'Fr': 223.0,
    'Ra': 226.0, 'Ac': 227.0, 'Th': 232.03811, 'Pa': 231.035882, 'U': 238.028913,
    'D': [-0.0368384615385, 0.928619230769],
    '1H': [-0.0368384615385, 0.928619230769]
}
