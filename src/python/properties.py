from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals
#
#   SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 
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
#	PROPERTIES
#
#	12/10/2009	--	initial coding				:	jc
#	11/24/2011	-- 	moved to seperate file	    :	jc
#	12/24/2015	-- 	refactored for release      :   jc
#	08/18/2016	-- 	doc strings                 :   jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Properties contains the classes that contain 
	atomic and molecular properties

'''
# OPEN	Need a consistent way to initialize "Selection/Atm" to avoid duplication


class Atomic(object):

    """ Class containing methods to define physical properties
        
        and atomic information used by other modules.

        http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some

        standard atomic weight is based on the natural istopic composition

        Examples
        ========

        >>> import sasmol.properties as properties
        >>> properties = properties.Atomic()
        >>> properties.amu()['FE']
        55.8452

        Note
        ----
    
        `self` parameter is not shown in the ``Parameters`` section in the documentation

    """ 

    def amu(self, **kwargs):
        """

        Parameters
        ----------
        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        standard_atomic_weights
            dictionary : containing standard atomic weights

        Examples
        -------

        >>> import sasmol.properties as properties
        >>> properties = properties.Atomic()
        >>> properties.amu()['FE']
        55.8452

        """

        mixed_standard_atomic_weight = {'H': 1.007947, 'He': 4.0026022,
                                        'Li': 6.9412, 'Be': 9.0121823, 'B': 10.8117, 'C': 12.01078,
                                        'N': 14.00672, 'O': 15.99943, 'F': 18.99840325, 'Ne': 20.17976,
                                        'Na': 22.9897702, 'Mg': 24.30506, 'Al': 26.9815382, 'Si': 28.08553,
                                        'P': 30.9737612, 'S': 32.0655, 'Cl': 35.4532, 'Ar': 39.9481,
                                        'K': 39.09831, 'Ca': 40.0784, 'Sc': 44.9559108, 'Ti': 47.8671,
                                        'V': 50.94151, 'Cr': 51.99616, 'Mn': 54.9380499, 'Fe': 55.8452,
                                        'Co': 58.9332009, 'Ni': 58.69342, 'Cu': 63.5463, 'Zn': 65.4094,
                                        'Ga': 69.7231, 'Ge': 72.641, 'As': 74.921602, 'Se': 78.963,
                                        'Br': 79.9041, 'Kr': 83.7982, 'Rb': 85.46783, 'Sr': 87.621,
                                        'Y': 88.905852, 'Zr': 91.2242, 'Nb': 92.906382, 'Mo': 95.942,
                                        'Tc': 98.0, 'Ru': 101.072, 'Rh': 102.905502, 'Pd': 106.421,
                                        'Ag': 107.86822, 'Cd': 112.4118, 'In': 114.8183, 'Sn': 118.7107,
                                        'Sb': 121.7601, 'Te': 127.603, 'I': 126.904473, 'Xe': 131.2936,
                                        'Cs': 132.905452, 'Ba': 137.3277, 'La': 138.90552, 'Ce': 140.1161,
                                        'Pr': 140.907652, 'Nd': 144.243, 'Pm': 145.0, 'Sm': 150.363,
                                        'Eu': 151.9641, 'Gd': 157.253, 'Tb': 158.925342, 'Dy': 162.5001,
                                        'Ho': 164.930322, 'Er': 167.2593, 'Tm': 168.93421, 'Yb': 173.043,
                                        'Lu': 174.9671, 'Hf': 174.9671, 'Ta': 180.94791, 'W': 183.841,
                                        'Re': 186.2071, 'Os': 190.233, 'Ir': 192.2173, 'Pt': 195.0782,
                                        'Au': 196.966552, 'Hg': 200.592, 'Tl': 204.38332, 'Pb': 207.21,
                                        'Bi': 208.980382, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0, 'Fr': 223.0,
                                        'Ra': 226.0, 'Ac': 227.0, 'Th': 232.03811, 'Pa': 231.035882,
                                        'U': 238.028913, 'D': 2.01410177804, '1H': 1.00782503214}
        #'U': 238.028913, 'D': 2.01410177804, '1H': 1.00782503214}

        standard_atomic_weight = {}
        for key, item in mixed_standard_atomic_weight.items():

            if 'keep_lower_case' in kwargs:
                standard_atomic_weight[key] = item
            else:
                standard_atomic_weight[key.upper()] = item

        return standard_atomic_weight



    def amino_acid_sld(self, **kwargs):
        """

        Parameters
        ----------
        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        residue_scattering
            dictionary : residue name, and sld information
            residue name : [vol Ang^3, eSL, SLprot Ang, SLdeut Ang, #exchngH]

        Examples
        -------

        >>> import sasmol.properties as properties
        >>> properties = properties.Atomic()
        >>> properties.amino_acid_sld()['PRO']
        [112.7, 52, 22.2207, 95.104, 0]

        Note
        ----

        Need to validate units here used in neutron reflectivity versus 
        small-angle scattering

        """
        
        residue_scattering = {
            'ALA': [88.6, 38, 20.1466,  61.7942, 1],
            'ARG': [173.4, 85, 56.9491, 129.8324, 6],
            'ASP': [111.1, 59, 42.149,  73.3816, 1],
            'ASN': [114.1, 60, 45.7009, 76.9366, 3],
            'CYS': [108.5, 54, 26.7345, 57.9702, 2],
            'GLU': [138.4, 67, 41.371,  93.372, 1],
            'GLN': [143.8, 68, 44.8675, 96.927, 3],
            'GLY': [60.1,  30, 20.98,   41.8038, 1],
            'HSD': [153.2, 72, 55.0709, 107.1304, 2],
            'HIS': [153.2, 72, 55.0709, 107.1304, 2],
            'HSE': [153.2, 72, 55.0709, 107.1304, 2],
            'HSP': [153.2, 72, 55.0709, 107.1304, 3],
            'ILE': [166.7, 62, 17.6464, 121.7654, 1],
            'LEU': [166.7, 62, 17.6464, 121.7654, 1],
            'LYS': [168.6, 71, 30.7473, 124.4544, 4],
            'MET': [162.9, 70, 21.3268, 104.622, 1],
            'PHE': [189.9, 78, 45.0734, 128.3686, 1],
            'PRO': [112.7, 52, 22.2207, 95.104, 0],
            'SER': [89.0,  46, 29.6925, 60.9282, 2],
            'THR': [116.1, 54, 25.1182, 87.5896, 2],
            'TRP': [227.8, 98, 67.7302, 151.0254, 2],
            'TYR': [193.6, 86, 54.6193, 127.5026, 2],
            'VAL': [140.0, 54, 18.4798, 101.775, 1]}

        return residue_scattering



    def set_average_vdw(self, **kwargs):

        """

        Parameters
        ----------
        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        self._atom_vdw updated in system object

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.set_average_vdw()
        >>> molecule.atom_vdw()[0]
        [-0.2, 1.85]

        Note
        ----

        Need to add citations.
    
        """

        element = self.element()

        van_der_waals = {'H': [-0.0368384615385, 0.928619230769],
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
                         '1H': [-0.0368384615385, 0.928619230769]}

        element_vdw = []

        for this_element in element:

            for key, item in van_der_waals.items():
                if this_element == key:
                    try:
                        if len(item) > 1:
                            element_vdw.append(item)
                    except:
                        element_vdw.append([None, None])

        self.setAtom_vdw(element_vdw)

        return

