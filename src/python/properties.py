from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals
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
# PROPERTIES
#
# 12/10/2009	--	initial coding				:	jc
# 11/24/2011	-- 	moved to seperate file	    :	jc
# 12/24/2015	-- 	refactored for release      :   jc
# 08/18/2016	-- 	doc strings                 :   jc
#
# 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
# *      **
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

    one_to_three_letter_protein_residue_dictionary = {
        'A': 'ALA',
        'R': 'ARG',
        'D': 'ASP',
        'N': 'ASN',
        'C': 'CYS',
        'E': 'GLU',
        'Q': 'GLN',
        'G': 'GLY',
        'H': 'HSE',
        'I': 'ILE',
        'L': 'LEU',
        'K': 'LYS',
        'M': 'MET',
        'F': 'PHE',
        'P': 'PRO',
        'S': 'SER',
        'T': 'THR',
        'W': 'TRP',
        'Y': 'TYR',
        'V': 'VAL',

    }

    one_to_three_letter_nucleic_residue_dictionary = {
        'G': 'GUA',
        'C': 'CYT',
        'A': 'ADE',
        'T': 'THY',
        'U': 'URA'
    }

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
        # 'U': 238.028913, 'D': 2.01410177804, '1H': 1.00782503214}

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

    def element_scattering_lengths(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        element_scattering_lengths
            dictionary : element name and scattering information

        """

        atomic_scattering = {
            'D': [2.014, 5.15, 0.282, 0.6671],
            'H': [1.008, 5.15, 0.282, -0.3741],
            'C': [12.01, 16.44, 1.692, 0.6646],
            'N': [14.01, 2.49, 1.974, 0.9360],
            'O': [16.00, 9.130, 2.256, 0.5803],
            'Na': [22.99, 4.45, 3.102, 0.3630],
            'Mg': [24.31, 1.56, 3.384, 0.5375],
            'K': [39.10, 11.01, 5.358, 0.3670],
            'Ca': [40.08, 4.19, 5.640, 0.4700],
            'Cl': [35.45, 24.84, 4.794, 0.9577],
            'Br': [79.90, 31.54, 9.870, 0.6795],
            'I': [126.9, 44.6, 14.946, 0.5280],
            'P': [30.97, 3.37, 4.230, 0.5130],
            'S': [32.07, 26.09, 4.512, 0.2847],
            'Fe': [55.85, 7.99, 7.332, 0.9450],
            'Co': [58.93, 7.99, 7.614, 0.2490],
            'Ni': [58.69, 8.18, 7.896, 1.0300],
            'Cu': [63.55, 8.78, 8.178, 0.7718],
            'Zn': [65.39, 9.85, 8.460, 0.5680],
        }

        return atomic_scattering

    def element_sl(self, **kwargs):
        return self.element_scattering_lengths(**kwargs)

    def nucleotide_scattering_lengths(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        nucleotide_scattering_lengths
            dictionary : nucleotide scattering information

        """

        residue_scattering = {
            'DA': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
            'DT': [304.0, 294.9, 44.2, 8.61, 9.65, 1, 12],
            'DG': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],
            'DC': [289.0, 280.3, 41.9, 8.68, 10.77, 2, 11],
            'U': [305.0, 296.9, 44.2, 9.28, 11.36, 2, 10],
            'A': [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
            'G': [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
            'C': [304.0, 295.9, 44.2, 8.68, 10.77, 2, 11],
        }

        return residue_scattering

    def nucleotide_sl(self, **kwargs):
        return self.nucleotide_scattering_lengths(**kwargs)

    def dna_scattering_lengths(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        dna_scattering_lengths
            dictionary : DNA scattering information

        """

        residue_scattering = {
            'DA': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
            'DT': [304.0, 294.9, 44.2, 8.61, 9.65, 1, 12],
            'DG': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],
            'DC': [289.0, 280.3, 41.9, 8.68, 10.77, 2, 11],
            'ADE': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
            'THY': [304.0, 294.9, 44.2, 8.61, 9.65, 1, 12],
            'GUA': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],
            'CYT': [289.0, 280.3, 41.9, 8.68, 10.77, 2, 11],
            'A': [313.0, 304.7, 45.3, 10.65, 12.73, 2, 11],
            'T': [304.0, 294.9, 44.2, 8.61, 9.65, 1, 12],
            'G': [329.0, 319.2, 47.6, 11.23, 14.35, 3, 11],
            'C': [289.0, 280.3, 41.9, 8.68, 10.77, 2, 11],
        }

        return residue_scattering

    def dna_sl(self, **kwargs):
        return self.dna_scattering_lengths(**kwargs)

    def rna_scattering_lengths(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        rna_scattering_lengths
            dictionary : RNA scattering information

        """

        residue_scattering = {
            'U': [305.0, 296.9, 44.2, 9.28, 11.36, 2, 10],
            'A': [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
            'G': [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
            'C': [304.0, 295.9, 44.2, 8.68, 10.77, 2, 11],
            'URA': [305.0, 296.9, 44.2, 9.28, 11.36, 2, 10],
            'ADE': [328.0, 319.2, 47.6, 10.65, 12.73, 2, 11],
            'GUA': [344.0, 334.9, 49.8, 11.23, 14.35, 3, 11],
            'CYT': [304.0, 295.9, 44.2, 8.68, 10.77, 2, 11],
        }

        return residue_scattering

    def rna_sl(self, **kwargs):
        return self.rna_scattering_lengths(**kwargs)

    def protein_scattering_lengths(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        protein_scattering_lengths
            dictionary : protein scattering information

        """

        residue_scattering = {
            'ALA': [71.1, 88.6, 10.7, 1.645, 2.686, 1, 5],
            'ARG': [156.2, 173.4, 24.0, 3.466, 9.714, 6, 13],
            'ASP': [115.1, 111.1, 16.7, 3.845, 4.886, 1, 4],
            'ASN': [114.1, 114.1, 16.9, 3.456, 6.580, 3, 6],
            'CYS': [103.1, 108.5, 15.2, 1.930, 4.013, 1, 5],
            'GLU': [129.1, 138.4, 18.9, 3.762, 4.803, 1, 6],
            'GLN': [128.1, 143.8, 19.2, 3.373, 6.497, 3, 8],
            'GLY': [57.1, 60.1, 8.5, 1.728, 2.769, 1, 3],
            'HSD': [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
            'HIS': [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
            'HSE': [137.2, 153.2, 20.2, 4.771, 6.521, 2, 7],
            'HSP': [138.2, 153.2, 20.2, 4.771, 6.521, 3, 7],
            'ILE': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
            'LEU': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
            'LYS': [128.2, 168.6, 20.1, 1.586, 5.752, 4, 13],
            'MET': [131.2, 162.9, 19.8, 1.763, 2.805, 1, 9],
            'PHE': [147.2, 189.9, 22.0, 4.139, 5.180, 1, 9],
            'PRO': [97.1, 112.7, 14.7, 2.227, 2.227, 0, 7],
            'SER': [87.1, 89.0, 13.0, 2.225, 4.308, 2, 5],
            'THR': [101.1, 116.1, 15.3, 2.142, 4.224, 2, 7],
            'TRP': [186.2, 227.8, 27.7, 6.035, 8.118, 2, 10],
            'TYR': [163.2, 193.6, 24.3, 4.719, 6.802, 2, 9],
            'VAL': [99.1, 140.0, 15.3, 1.478, 2.520, 1, 9],
            'A': [71.1, 88.6, 10.7, 1.645, 2.686, 1, 5],
            'R': [156.2, 173.4, 24.0, 3.466, 9.714, 6, 13],
            'D': [115.1, 111.1, 16.7, 3.845, 4.886, 1, 4],
            'N': [114.1, 114.1, 16.9, 3.456, 6.580, 3, 6],
            'C': [103.1, 108.5, 15.2, 1.930, 4.013, 1, 5],
            'E': [129.1, 138.4, 18.9, 3.762, 4.803, 1, 6],
            'Q': [128.1, 143.8, 19.2, 3.373, 6.497, 3, 8],
            'G': [57.1, 60.1, 8.5, 1.728, 2.769, 1, 3],
            'H': [138.2, 153.2, 20.2, 4.771, 6.521, 3, 7],
            'I': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
            'L': [113.2, 166.7, 17.5, 1.395, 2.437, 1, 11],
            'K': [128.2, 168.6, 20.1, 1.586, 5.752, 4, 13],
            'M': [131.2, 162.9, 19.8, 1.763, 2.805, 1, 9],
            'F': [147.2, 189.9, 22.0, 4.139, 5.180, 1, 9],
            'P': [97.1, 112.7, 14.7, 2.227, 2.227, 0, 7],
            'S': [87.1, 89.0, 13.0, 2.225, 4.308, 2, 5],
            'T': [101.1, 116.1, 15.3, 2.142, 4.224, 2, 7],
            'W': [186.2, 227.8, 27.7, 6.035, 8.118, 2, 10],
            'Y': [163.2, 193.6, 24.3, 4.719, 6.802, 2, 9],
            'V': [99.1, 140.0, 15.3, 1.478, 2.520, 1, 9],
        }

        return residue_scattering

    def protein_sl(self, **kwargs):
        return self.protein_scattering_lengths(**kwargs)

    def van_der_waals_radii(self, **kwargs):
        """

        Parameters
        ----------
        kwargs
            optional future arguments

        Returns
        -------
        van_der_waals_radii
            dictionary : element name and van der Waals radius

        """

        non_bonded_vdw = {
            'H': 1.20, 'He': 1.40,
            'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.70,
            'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54,
            'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.1,
            'P': 1.8, 'S': 1.8, 'Cl': 1.75, 'Ar': 1.88,
            'K': 2.75, 'Ca': 2.31, 'Sc': 2.11, 'Ti': 2.15,
            'V': 2.05, 'Cr': 2.05, 'Mn': 2.05, 'Fe': 2.05,
            'Co': 2.0, 'Ni': 1.63, 'Cu': 1.4, 'Zn': 1.39,
            'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.9,
            'Br': 1.85, 'Kr': 2.02, 'Rb': 3.03, 'Sr': 2.49,
            'Y': 2.4, 'Zr': 2.3, 'Nb': 2.15, 'Mo': 2.1,
            'Tc': 2.05, 'Ru': 2.05, 'Rh': 2.0, 'Pd': 1.63,
            'Ag': 1.72, 'Cd': 1.58, 'In': 1.93, 'Sn': 2.17,
            'Sb': 2.06, 'Te': 2.06, 'I': 1.98, 'Xe': 2.16,
            'Cs': 3.43, 'Ba': 2.68, 'La': 2.5, 'Ce': None,
            'Pr': None, 'Nd': None, 'Pm': None, 'Sm': None,
            'Eu': None, 'Gd': None, 'Tb': None, 'Dy': None,
            'Ho': None, 'Er': None, 'Tm': None, 'Yb': None,
            'Lu': None, 'Hf': 2.25, 'Ta': 2.2, 'W': 2.1,
            'Re': 2.05, 'Os': 2.0, 'Ir': 2.0, 'Pt': 1.75,
            'Au': 1.66, 'Hg': 1.55, 'Tl': 1.96, 'Pb': 2.02,
            'Bi': 2.07, 'Po': 1.97, 'At': 2.02, 'Rn': 2.20,
            'Fr': 3.48, 'Ra': 2.83, 'Ac': None, 'Th': 2.4,
            'Pa': None, 'U': 1.86, 'D': 1.2, '1H': 1.2,
        }

        vdw = {}
        for key, item in non_bonded_vdw.items():
            if 'keep_lower_case' in kwargs:
                vdw[key] = item
            else:
                vdw[key.upper()] = item

        return vdw

    def van_der_Waals_radii(self, **kwargs):
        return self.van_der_waals_radii(**kwargs)

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
