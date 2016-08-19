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

    def charmm_names(self, **kwargs):
        """

        Parameters
        ----------
        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        charmm_names
            dictionary : atom_name and list of atom names for
            hydrogen, carbon, nitrogen, oxygen, sulfur, phosphorus, other

        Examples
        -------

        >>> import sasmol.properties as properties
        >>> properties = properties.Atomic()
        >>> properties.charmm_names()['sulfur']
        ['SG', 'SD', '1SG', '2SG']

        """

        charmm_names = {}

        charmm_names['hydrogen'] = ['HN', 'HA', 'HB1', 'HB2', 'HB3', 'HG1', 'HG2', 'HD1', 'HD2', 'HE', 'HH11', 'HH12', 'HH21', 'HH22', 'HD21', 'HD22', 'HE21', 'HE22', 'HA1', 'HA2', 'HE1', 'HE2', 'HB', 'HG21', 'HG22', 'HG23', 'HG11', 'HG12', 'HD3', 'HG', 'HD11', 'HD12', 'HD13', 'HD23', 'HZ1', 'HZ2', 'HZ3', 'HE3', 'HZ', 'HH2', 'HH', 'HG13', 'H1', 'H2', 'HC', 'HD', 'HMA1', 'HMA2', 'HMA3', 'HAA1', 'HAA2', 'HBA1', 'HBA2', 'HMB1', 'HMB2', 'HMB3', 'HAB', 'HBB1', 'HBB2', 'HMC1', 'HMC2', 'HMC3', 'HAC', 'HBC1', 'HBC2', 'HMD1', 'HMD2', 'HMD3', 'HAD1', 'HAD2', 'HBD1', 'HBD2', 'HT1', 'HT2', 'HT3', 'HN1', 'HN2', 'HY1', 'HY2',
                    'HY3', 'HNT', "H5'", "H5''", "H4'", "H1'", 'H21', 'H22', 'H8', "H2''", "H2'", "H3'", 'H61', 'H62', 'H6', 'H5', 'H41', 'H42', 'H3', 'H51', 'H52', 'H53', 'H11', 'H12', 'H13', 'H23', 'H9', 'H5T', "H53'", 'H5T1', 'H5T2', 'H5T3', 'H3T', 'H3T1', 'H3T2', 'H3T3', 'H91', 'H92', 'H93', 'H9B1', 'H9B2', 'H9B3', 'H5M1', 'H5M2', 'H5M3', 'H4', 'H15', 'H16', 'H17', "H11'", "H12'", "H21'", "H22'", "H31'", "H32'", "H41'", "H42'", "H51'", "H52'", "H1''", "H4''", "H3''", "1H3'", "1H3''", "1H2''", "2H5'", "2H5''", "AH4'", "AH1'", 'AH8', 'AH2', "AH2'", 'AH61', 'AH62', 'AH2T', "AH3'", 'AH3T', "AH5'", 'AH5S']

        charmm_names['carbon'] = ["CA", "CB", "C", "CG", "CD", "CZ", "CE1", "CD2", "CG2", "CG1", "CD1", "CE", "CE2", "CE3", "CZ3", "CZ2", "CH2", "C1A", "C2A", "C3A", "C4A", "C1B", 'C2B', 'C3B', 'C4B', 'C1C', 'C2C', 'C3C', 'C4C', 'C1D', 'C2D', 'C3D', 'C4D', 'CHA', 'CHB', 'CHC', 'CHD', 'CMA', 'CAA', 'CBA', 'CGA', 'CMB', 'CAB', 'CBB', 'CMC',
                  'CAC', 'CBC', 'CMD', 'CAD', 'CBD', 'CGD', 'CAY', 'CY', 'CT', 'CAT', "C5'", "C4'", "C1'", 'C4', 'C2', 'C6', 'C5', 'C8', "C2'", "C3'", 'C5M', 'C1', 'C5T', 'C3T', 'C9', 'C9B', 'C3', 'C7', 'C10', 'C12', '1CB', '2CB', "1C3'", "1C2'", "2C5'", "AC4'", "AC1'", 'AC5', 'AC8', 'AC2', 'AC4', 'AC6', "AC2'", "AC3'", "AC5'"]

        charmm_names['nitrogen'] = ['N', 'NH1', 'NH2', 'ND2', 'NE2', 'ND1', 'NZ', 'NE1', 'NA', 'NB', 'NC', 'ND', 'NT', 'N9', 'N2', 'N3', 'N1', 'N7', 'N6', 'N4', 'N14', 'NP', 'NO1', 'NO2', "NO5'", "NC5'", 'NH5S', "NH5'", "NC2'", "NH2'", "NO2'",
                    'NH2T', "NC3'", "NH3'", "NO3'", 'NH3T', "NC1'", "NH1'", "NC4'", "NH4'", "NO4'", 'NN1', 'NC6', 'NH6', 'NC5', 'NH5', 'NC4', 'NH4', 'NC3', 'NC2', 'NC7', 'NO7', 'NN7', 'NH71', 'NH72', 'NH42', 'AN7', 'AN9', 'AN1', 'AN3', 'AN6']

        charmm_names['oxygen'] = ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'OH2', 'O1A', 'O2A', 'O1D', 'O2D', 'OY', 'OT1', 'OT2', 'O1', 'O2', 'O3', 'O4', 'O1P', 'O2P', "O5'", "O4'", 'O6', "O2'", "O3'", 'O5T', 'O1P3', 'O2P3', 'O3T', 'O3P3', 'O13',
                  'O11', 'O12', 'O14', 'O22', 'O23', 'O24', 'O3A', 'O1B', 'O2B', 'O3B', 'O1G', 'O2G', 'O3G', 'O21', 'O3P', 'O31P', 'O32P', 'O33P', "1O2'", '2O1P', '2O2P', "2O5'", "AO4'", "AO2'", "AO3'", 'AO1', 'AO2', "AO5'", 'AO1P', 'AO2P', 'AO2T', 'OM']

        charmm_names['sulfur'] = ['SG', 'SD', '1SG', '2SG']

        charmm_names['phosphorus'] = ['P1', 'P', 'P3', 'P2',
                      'PA', 'PB', 'PG', '2P', 'AP', 'AP2']

        charmm_names['other'] = ['AS', 'AG', 'AL', 'AR', 'AC', 'AT', 'AU', 'BI', 'BE', 'B', 'BR', 'BA', 'CR', 'CS', 'CU', 'CO', 'D', 'DY', 'EU', 'ER', 'F', 'FE', 'FR', 'GA', 'GE', 'GD', 'HO', 'HF', 'IN', 'I', 'IR', 'K', 'KR', 'LI', 'LA', 'LU', 'MG', 'MN', 'MO', 'NE',
                 'NI', 'OS', 'PD', 'PR', 'PM', 'PT', 'PO', 'RB', 'RU', 'RH', 'RE', 'RN', 'RA', 'SI', 'SC', 'SE', 'SR', 'SN', 'SB', 'SM', 'TI', 'TC', 'TE', 'TB', 'TA', 'TL', 'TH', 'U', 'V', 'W', 'XE', 'Y', 'YB', 'ZN', 'ZR', 'SS', 'CAL', 'DUM', 'POT', 'CES', 'CLA']

# note that SM is C-S-S-C not Sm (the element)
# note that SM is C-S-S-C not Sm (the element)

        return charmm_names

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

    def create_fasta(self, **kwargs):
        """
        Method to make a fasta file compatible lists of residue names

        http://en.wikipedia.org/wiki/FASTA_format

        Parameters
        ----------
        kwargs 
            optional future arguments
                                                                                     
        Returns
        -------
        self._fasta
            sets the fasta attribute in a system object 

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule("hiv1_gag.pdb")
        >>> molecule.create_fasta()
        >>> molecule.fasta()[:5]
        ['G', 'A', 'R', 'A', 'S']
        >>> m.create_fasta(fasta_format=True)
        >>> print molecule.fasta()[:5]
        >
        GAR
        >>> m.create_fasta(fasta_format=True,width='60')
        >>> print molecule.fasta()
        >
        GARASVLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQIL
        GQLQPSLQTGSEELRSLYNTIAVLYCVHQRIDVKDTKEALDKIEEEQNKSKKKAQQAAAD
        TGNNSQVSQNYPIVQNLQGQMVHQAISPRTLNAWVKVVEEKAFSPEVIPMFSALSEGATP
        QDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRVHPVHAGPIAPGQMREPRGSDIAGTTS
        TLQEQIGWMTNNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFY
        KTLRAEQASQEVKNAATETLLVQNANPDCKTILKALGPAATLEEMMTACQGVGGPGHKAR
        VIAEAMSQVTNSATIMMQKGNFRNQRKTVKCFNCGKEGHIAKNCRAPRKKGCWKCGKEGH
        QMKDCTERQAN

        >>> molecule.create_fasta(fasta_format=True,width='60',name='aar')
        >>> print molecule.fasta()[:90]
        >aar
        GARASVLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQIL
        GQLQPSLQTGSEELRSLYNTIAVL

        
        Note
        ----
        Other use types below

        molecule.create_fasta(fasta_format=True,exclude_hetatm=True,by_chain=True)
        print molecule.fasta()

        print '>>> testing by_chain with HETATM (t.py): '

        molecule.create_fasta(fasta_format=True,by_chain=True)
        print molecule.fasta()

        print '>>> testing by_segname (t.py): '

        molecule.create_fasta(fasta_format=True,exclude_hetatm=True,by_segname=True)
        print molecule.fasta()

        Note that this creates a simple string that is associated with the molecule (self)
        and it will return without assigning a string if a non-standard three letter code
        is supplied.

        """

        fasta_format = False
        by_chain = False
        by_segname = False
        single_chain = False
        single_segname = False
        exclude_hetatm = False

        if 'fasta_format' in kwargs:
            fasta_format = True

        if 'name' in kwargs:
            name = kwargs['name']
            header = '>' + name
        else:
            header = '>'

        if 'width' in kwargs:
            width = kwargs['width']
        else:
            width = '80'

        if 'exclude_hetatm' in kwargs:
            exclude_hetatm = True

        if 'by_chain' in kwargs:
            by_chain = True

        elif 'by_segname' in kwargs:
            by_segname = True

        one_resname = []

        resname = self.resname()
        atom = self.atom()
        chain = self.chain()
        segname = self.segname()

        residue_dictionary = {
            'ALA': 'A',
            'ARG': 'R',
            'ASP': 'D',
            'ASN': 'N',
            'CYS': 'C',
            'GLU': 'E',
            'GLN': 'Q',
            'GLY': 'G',
            'HSD': 'H',
            'HIS': 'H',
            'HSE': 'H',
            'HSP': 'H',
            'ILE': 'I',
            'LEU': 'L',
            'LYS': 'K',
            'MET': 'M',
            'PHE': 'F',
            'PRO': 'P',
            'SER': 'S',
            'THR': 'T',
            'TRP': 'W',
            'TYR': 'Y',
            'VAL': 'V',
            'GUA': 'G',
            'CYT': 'C',
            'ADE': 'A',
            'THY': 'T',
            'URA': 'U'
        }

        for i in xrange(len(resname)):
            this_resname = resname[i]
            if this_resname in residue_dictionary:
                one_resname.append(residue_dictionary[this_resname])
            elif (atom[i] == 'HETATM'):
                # print 'skipping non-standard resname in HETATM record:
                # ',this_resname
                one_resname.append("X")
            else:
                print('non standard resname: ', this_resname)
                print('unable to make Fasta conversion')
                return

        self.setOne_letter_resname(one_resname)

        resid = self.resid()

        last_resid = None
        last_chain = None
        last_segname = None
        fasta = []
        local_fasta = []
        number_of_chains = 0
        number_of_segnames = 0
        chain_name = []
        segname_name = []
        first = True
        for i in xrange(len(resid)):
            this_resid = resid[i]

            if by_chain:
                this_chain = chain[i]
                if this_chain != last_chain:
                    number_of_chains += 1
                    chain_name.append(this_chain)
                    if first:
                        first = False
                    else:
                        fasta.append(local_fasta)

                    last_chain = this_chain
                    local_fasta = []

            if by_segname:
                this_segname = segname[i]
                if this_segname != last_segname:
                    number_of_segnames += 1
                    segname_name.append(this_segname)
                    if first:
                        first = False
                    else:
                        fasta.append(local_fasta)

                    last_segname = this_segname
                    local_fasta = []

            this_resname = self.one_letter_resname()[i]
            if this_resid != last_resid:
                local_fasta.append(this_resname)
                last_resid = this_resid

        if by_chain or by_segname:
            fasta.append(local_fasta)

        elif not by_chain and not by_segname:
            fasta = local_fasta
            number_of_chains = 1
            chain_name = chain[0]

        final_fasta = ''

        if by_segname:
            number_of_chains = number_of_segnames

        for i in xrange(number_of_chains):
            saveme = False
            if fasta_format:
                from textwrap import TextWrapper
                wrapper = TextWrapper(width=int(width))
                if by_chain or by_segname:
                    if exclude_hetatm:
                        while "X" in fasta[i]:
                            fasta[i].remove("X")

                    joined_fasta = ''.join(fasta[i])
                else:
                    if exclude_hetatm:
                        while "X" in fasta:
                            fasta.remove("X")
                    joined_fasta = ''.join(fasta)

                formatted_fasta = "\n".join(wrapper.wrap(joined_fasta))

                if(len(formatted_fasta.strip()) > 0):
                    saveme = True

                if(len(header) > 1):
                    if by_chain:
                        formatted_fasta = header + ' chain:' + \
                            chain_name[i] + '\n' + formatted_fasta
                    elif by_segname:
                        formatted_fasta = header + ' segname:' + \
                            segname_name[i] + '\n' + formatted_fasta
                    else:
                        formatted_fasta = header + '\n' + formatted_fasta

                else:
                    if by_chain:
                        formatted_fasta = header + 'chain:' + \
                            chain_name[i] + '\n' + formatted_fasta
                    elif by_segname:
                        formatted_fasta = header + 'segname:' + \
                            segname_name[i] + '\n' + formatted_fasta
                    else:
                        formatted_fasta = header + '\n' + formatted_fasta

                if saveme:
                    final_fasta += formatted_fasta + '\n'
                else:
                    final_fasta += '\n'

        if fasta_format:
            self.setFasta(final_fasta)
        else:
            self.setFasta(fasta)

        return

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

