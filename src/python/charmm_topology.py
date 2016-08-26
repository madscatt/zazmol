from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

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
#	TOPOLOGY
#
#	1/26/2012	--	initial coding				: jc
#	12/26/2015	--	refactored for release      : jc
#	08/18/2016	--	documentation               : jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **

import os
import numpy
import copy
import sasmol.config as config

class CharmmTopology(object):

    '''
    This class contains charmm topology information used other modules.


    The output topology dictionary looks like the following examples:

    Examples
    ========

    >>> import pprint
    >>> import sasmol.topology as topology
    >>> test = topology.CharmmTopology()
    >>> test.read_charmm_topology()
    >>> pprint.pprint(test.topology_info['NTER'],width=100)
    {'ATOM': [['N', 'NH3', '-0.30'],
          ['HT1', 'HC', '0.33'],
          ['HT2', 'HC', '0.33'],
          ['HT3', 'HC', '0.33'],
          ['CA', 'CT1', '0.21'],
          ['HA', 'HB', '0.10']],
    'BOND': [['HT1', 'N'], ['HT2', 'N'], ['HT3', 'N']],
    'DELE': [['ATOM', 'HN']],
    'DONO': ['HT1', 'N', 'HT2', 'N', 'HT3', 'N'],
    'IC': [['HT1', 'N', 'CA', 'C', '0.0000', '0.0000', '180.0000', '0.0000', '0.0000'],
            ['HT2', 'CA', '*N', 'HT1', '0.0000', '0.0000', '120.0000', '0.0000', '0.0000'],
            ['HT3', 'CA', '*N', 'HT2', '0.0000', '0.0000', '120.0000', '0.0000', '0.0000']],
    'TOTAL_CHARGE': '1.00'}
    
    >>> import pprint
    >>> import sasmol.topology as topology
    >>> test = topology.CharmmTopology()
    >>> test.read_charmm_topology()
    >>> pprint.pprint(test.topology_info['ALA'],width=100)
    {'ACCE': ['O', 'C'],
    'ATOM': [['N', 'NH1', '-0.47'],
            ['HN', 'H', '0.31'],
            ['CA', 'CT1', '0.07'],
            ['HA', 'HB', '0.09'],
            ['CB', 'CT3', '-0.27'],
            ['HB1', 'HA', '0.09'],
            ['HB2', 'HA', '0.09'],
            ['HB3', 'HA', '0.09'],
            ['C', 'C', '0.51'],
            ['O', 'O', '-0.51']],
    'BOND': [['CB', 'CA'],
            ['N', 'HN'],
            ['N', 'CA'],
            ['C', 'CA'],
            ['C', '+N'],
            ['CA', 'HA'],
            ['CB', 'HB1'],
            ['CB', 'HB2'],
            ['CB', 'HB3']],
    'DONO': ['HN', 'N'],
    'DOUB': [['O', 'C']],
    'IC': [['-C', 'CA', '*N', 'HN', '1.3551', '126.4900', '180.0000', '115.4200', '0.9996'],
            ['-C', 'N', 'CA', 'C', '1.3551', '126.4900', '180.0000', '114.4400', '1.5390'],
            ['N', 'CA', 'C', '+N', '1.4592', '114.4400', '180.0000', '116.8400', '1.3558'],
            ['+N', 'CA', '*C', 'O', '1.3558', '116.8400', '180.0000', '122.5200', '1.2297'],
            ['CA', 'C', '+N', '+CA', '1.5390', '116.8400', '180.0000', '126.7700', '1.4613'],
            ['N', 'C', '*CA', 'CB', '1.4592', '114.4400', '123.2300', '111.0900', '1.5461'],
            ['N', 'C', '*CA', 'HA', '1.4592', '114.4400', '-120.4500', '106.3900', '1.0840'],
            ['C', 'CA', 'CB', 'HB1', '1.5390', '111.0900', '177.2500', '109.6000', '1.1109'],
            ['HB1', 'CA', '*CB', 'HB2', '1.1109', '109.6000', '119.1300', '111.0500', '1.1119'],
            ['HB1', 'CA', '*CB', 'HB3', '1.1109', '109.6000', '-119.5800', '111.6100', '1.1114']],
    'IMPR': [['N', '-C', 'CA', 'HN'], ['C', 'CA', '+N', 'O']],


        Note
        ----
    
        `self` parameter is not shown in the ``Parameters`` section in the documentation


    '''

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


    def add(self, dictionary, key, value, **kwargs):

        '''
        This method is used by the self.read_topology method.

        It will check whether the key is in the dictionary:
        if yes, append value to dictionary[key]
        if no, initialize dictionary[key] as [value]

        Parameters
        ----------
        dictionary 
            dictionary : containing existing key, value pairs

        key 
            dictionary key : key to query and perhaps intialize

        value 
            dictionary value : value add

        kwargs 
            optional future arguments


        Returns
        -------
        updated dictionary


        Examples
        -------

        >>> import sasmol.topology as topology
        >>> d = {}
        >>> d['ALA'] = 'RESI ALA'
        >>> test = topology.CharmmTopology()
        >>> test.add(d, 'GLU', 'GLU')
        >>> d
        {'GLU': ['GLU'], 'ALA': 'RESI ALA'} 
        >>> test.add(d, 'ASP', 'RESI ASP')
        >>> d
        {'ASP': ['RESI ASP'], 'GLU': ['GLU'], 'ALA': 'RESI ALA'}
       
        ''' 
   
        if key not in dictionary:
            dictionary[key] = [value]
        else:
            dictionary[key].append(value)


    def read_charmm_topology(self, topology_file_path='', topology_file_name='top_all27_prot_na.inp', **kwargs):
        '''
        Read and parse the charmm topology file
        A comprehensive dictionary (topology_info) will be built to store all the topology information
        The strategy for parsing topology file is to split words in each line

        Parameters
        ----------
        topology_file_path 

            string : default = ''

        topology_file_name 

            string : default = 'top_all27_prot_na.inp'

        kwargs 
            optional future arguments

        Examples
        -------

        >>> import pprint
        >>> import sasmol.topology as topology
        >>> test = topology.CharmmTopology()
        >>> test.read_charmm_topology()
        >>> pprint.pprint(test.topology_info['NTER'],width=100)
        {'ATOM': [['N', 'NH3', '-0.30'],
            ['HT1', 'HC', '0.33'],
            ['HT2', 'HC', '0.33'],
            ['HT3', 'HC', '0.33'],
            ['CA', 'CT1', '0.21'],
            ['HA', 'HB', '0.10']],
        'BOND': [['HT1', 'N'], ['HT2', 'N'], ['HT3', 'N']],
        'DELE': [['ATOM', 'HN']],
        'DONO': ['HT1', 'N', 'HT2', 'N', 'HT3', 'N'],
        'IC': [['HT1', 'N', 'CA', 'C', '0.0000', '0.0000', '180.0000', '0.0000', '0.0000'],
                ['HT2', 'CA', '*N', 'HT1', '0.0000', '0.0000', '120.0000', '0.0000', '0.0000'],
                ['HT3', 'CA', '*N', 'HT2', '0.0000', '0.0000', '120.0000', '0.0000', '0.0000']],
        'TOTAL_CHARGE': '1.00'}
    
        Returns
        -------
        error
            list : with error code

        Note
        ----
    
        error code needs to be addressed and handled globally


        '''
        error = []
        self.topology_info = {}
        lines = open(os.path.join(topology_file_path,
                                  topology_file_name), 'r').readlines()
        for line in lines:
            line = line.strip()
            if (len(line) > 0 and line[0] != '!'):
                words = line.split()
                if words[0][0:4] == 'MASS':
                    self.add(self.topology_info, 'MASS', words[1:4])
                #
                elif words[0][0:4] == 'DECL':
                    for ind in range(1, 2):
                        self.add(self.topology_info, 'DECL', words[ind])
                #
                elif words[0][0:4] == 'DEFA':  # Need further parsing
                    for ind in range(1, 5):
                        self.add(self.topology_info, 'DEFA', words[ind])
                #
                elif words[0][0:4] == 'AUTO':  # Need further parsing
                    for ind in range(1, 3):
                        self.add(self.topology_info, 'AUTO', words[ind])
                #
                elif words[0][0:4] == 'RESI' or words[0][0:4] == 'PRES':
                    cur_res = words[1]
                    self.topology_info[cur_res] = {}
                    self.topology_info[cur_res]['TOTAL_CHARGE'] = words[2]
                #
                else:
                    # elif 'cur_res' in locals():
                    #
                    if words[0] == 'ATOM':
                        self.add(self.topology_info[
                                 cur_res], 'ATOM', words[1:4])
                    #
                    elif words[0][0:4] == 'BOND':
                        ind = 1
                        while (ind + 2 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the BOND line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'BOND', words[ind:ind + 2])
                            ind += 2
                    #
                    elif words[0][0:4] == 'DOUB':
                        ind = 1
                        while (ind + 2 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the DOUBLE line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'DOUB', words[ind:ind + 2])
                            ind += 2
                    #
                    elif words[0][0:4] == 'IMPR':
                        ind = 1
                        while (ind + 4 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])
                                    and (words[ind + 3][0].isalnum() or words[ind + 3][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the IMPR line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'IMPR', words[ind:ind + 4])
                            ind += 4
                    #
                    elif words[0][0:4] == 'CMAP':
                        ind = 1
                        while (ind + 4 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])
                                    and (words[ind + 3][0].isalnum() or words[ind + 3][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the CMAP line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'CMAP', words[ind:ind + 4])
                            ind += 4
                    #
                    elif words[0][0:4] == 'ANGL':
                        ind = 1
                        while (ind + 3 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the ANGLE line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'ANGL', words[ind:ind + 3])
                            ind += 3
                    #
                    elif words[0][0:4] == 'THET':
                        ind = 1
                        while (ind + 3 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the THET line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'THET', words[ind:ind + 3])
                            ind += 3
                    #
                    elif words[0][0:4] == 'DIHE':
                        ind = 1
                        while (ind + 4 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])
                                    and (words[ind + 1][0].isalnum() or words[ind + 1][0] in ['+', '-'])
                                    and (words[ind + 2][0].isalnum() or words[ind + 2][0] in ['+', '-'])
                                    and (words[ind + 3][0].isalnum() or words[ind + 3][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the DIHE line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'DIHE', words[ind:ind + 4])
                            ind += 4
                    #
                    elif words[0][0:4] == 'DONO':  # Need to check
                        ind = 1
                        while (ind + 1 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the DONOR line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'DONO', words[ind])
                            ind += 1
                    #
                    elif words[0][0:4] == 'ACCE':  # Need to check
                        ind = 1
                        while (ind + 1 <= len(words) and words[ind][0] != '!'):
                            if not ((words[ind][0].isalnum() or words[ind][0] in ['+', '-'])):
                                error.append(
                                    '\nSomething wrong with the ACCEPTOR line: \n' + line + '\nof residue: ' + cur_res + '\n')
                                return error
                            self.add(self.topology_info[
                                     cur_res], 'ACCE', words[ind])
                            ind += 1
                    #
                    elif words[0][0:4] == 'IC':
                        self.add(self.topology_info[
                                 cur_res], 'IC', words[1:10])
                    #
                    elif words[0][0:4] == 'BUIL':
                        self.add(self.topology_info[
                                 cur_res], 'BUIL', words[1:10])
                    #
                    elif words[0][0:4] == 'PATC':  # Need further parsing
                        self.add(self.topology_info[
                                 cur_res], 'PATC', words[1:])
                    #
                    elif words[0][0:4] == 'DELE':
                        if 'DELE' not in self.topology_info[cur_res]:
                            self.topology_info[cur_res]['DELE'] = {}
                        prop = words[1][0:4]
                        if prop == 'ATOM':
                            # Only 1 word  is parsed after 'DELE ATOM', which
                            # is ok for "top_all27_prot_na.inp"
                            self.add(self.topology_info[cur_res][
                                     'DELE'], 'ATOM', words[2])
                        elif prop == 'ANGL':
                            # Need further parsing
                            self.add(self.topology_info[cur_res][
                                     'DELE'], 'ANGL', words[2:])
                        else:
                            # print 'WARNING!\n'+line+' NOT PARSED!\n' #
                            # Current "top_all27_prot_na.inp" file only has
                            # DELE ATOM and DELE ANGL, and therefore other
                            # situations are not parsed but warned
                            pass
        return error

    def patch_charmm_residue_atoms(self, residue, patch, **kwargs):

        '''
        Applies 'patch' to 'residue' based on the charmm topology 
        definition to add / remove atoms as needed.  Only ATOM and DELE
        records are accommodated
        
        Parameters
        ----------
        residue 
            string : residue name

        patch
            string : patch name

        kwargs 
            optional future arguments

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.check_charmm_atomic_order_reorganize()
        >>> molecule.patch_charmm_residue_atoms('GLY','GLYP') 
        
        Returns
        -------
        updated dictionary entry for residue in question

        Note
        ____

        Error list is defined and used, but not returned.  Needs to be handled globally.

            
        '''
        
        # OPEN: ONLY ATOMS ARE PATCHED (NO BOND, ETC)
        error = []
        #
        # Make sure residue and patch are the right type
        #
        if residue not in self.topology_info:
            error.append('Residue ' + residue +
                         ' not found in the topology information list during patching!')
            return error
        elif patch not in self.topology_info:
            error.append(
                'Patch ' + residue + ' not found in the topology information list during patching!')
            return error
        #
        # Delete atoms from 'DELE' in patch list
        #
        tmp_dict = copy.deepcopy(self.topology_info[residue])
        # print 'ZHL ',patch,self.topology_info[patch]
        if 'DELE' in self.topology_info[patch]['DELE']:
            for atom_delete in self.topology_info[patch]['DELE']['ATOM']:
                for atom_residue in tmp_dict['ATOM']:
                    if atom_delete == atom_residue[0]:
                        tmp_dict['ATOM'].remove(atom_residue)
        #
        # Delete atoms
        #
        num_added = 0
        for atom_patch in self.topology_info[patch]['ATOM']:
            for atom_residue in tmp_dict['ATOM']:
                if atom_residue[0] == atom_patch[0]:
                    tmp_dict['ATOM'].remove(atom_residue)
            if patch in ['NTER', 'GLYP', 'PROP']:
                tmp_dict['ATOM'].insert(num_added, atom_patch)
            if patch == 'CTER':
                tmp_dict['ATOM'].append(atom_patch)
            num_added += 1
        #
        # Add the patched residue to the topology dictionary as a new key of eg 'ALA_NTER'
        #
        self.topology_info[residue + '_' + patch] = tmp_dict
        self.charmm_residue_atoms[
            residue + '_' + patch] = numpy.array(tmp_dict['ATOM'])[:, 0].tolist()

    def setup_cys_patch_atoms_simple(self, **kwargs):
        '''
        A way to set up cys patch atoms
        simply remove HG1 atom from the atom list

        Parameters
        ----------
        
        kwargs 
            optional future arguments

        Examples
        -------

        >>> import sasmol.topology as topology
        >>> import pprint
        >>> test = topology.CharmmTopology()
        >>> test.read_charmm_topology()
        >>> pprint.pprint(self.topology_info['CYS']['ATOM'])
        [['N', 'NH1', '-0.47'],
        ['HN', 'H', '0.31'],
        ['CA', 'CT1', '0.07'],
        ['HA', 'HB', '0.09'],
        ['CB', 'CT2', '-0.11'],
        ['HB1', 'HA', '0.09'],
        ['HB2', 'HA', '0.09'],
        ['SG', 'S', '-0.23'],
        ['HG1', 'HS', '0.16'],
        ['C', 'C', '0.51'],
        ['O', 'O', '-0.51']]
        >>> pprint.pprint(test.setup_cys_patch_atoms_simple(), width=80)
        [['N', 'NH1', '-0.47'],
        ['HN', 'H', '0.31'],
        ['CA', 'CT1', '0.07'],
        ['HA', 'HB', '0.09'],
        ['CB', 'CT2', '-0.11'],
        ['HB1', 'HA', '0.09'],
        ['HB2', 'HA', '0.09'],
        ['SG', 'S', '-0.23'],
        ['C', 'C', '0.51'],
        ['O', 'O', '-0.51']] 

        
        Returns
        -------
        updated dictionary entry for residue in question

        '''

        atoms = self.topology_info['CYS']['ATOM']
        for atom in atoms:
            if atom[0] == 'HG1':
                atoms.remove(atom)
        # print atoms
        return atoms

    def setup_charmm_residue_atoms(self, **kwargs):
        '''

        Build the atom list of all the residues in the charmm topology file

        Parameters
        ----------
        
        kwargs 
            optional future arguments

        Examples
        -------

        >>> import sasmol.topology as topology
        >>> test = topology.CharmmTopology()
        >>> test.read_charmm_topology()
        >>> test.setup_charmm_residue_atoms() 

        Returns
        -------
        self.charmm_residue_atoms : new dictionary entry for residue in question

        '''

        self.charmm_residue_atoms = {}
        for (key, value) in zip(self.topology_info.keys(), self.topology_info.values()):
            if type(value) is dict:
                if 'ATOM' in value:
                    atoms = numpy.array(value['ATOM'])[:, 0].tolist()
                    self.charmm_residue_atoms[key] = atoms
        # OPEN Setup CYS patch atoms in a simple way
        self.charmm_residue_atoms['DISU'] = self.setup_cys_patch_atoms_simple()

    def compare_list_ignore_order(self, l1, l2, **kwargs):
        '''
        Compare two lists while ignoring order

        Parameters
        ----------
        l1
            list : list 1 

        l2
            list : list 2

        kwargs 
            optional future arguments

        Examples
        -------

        >>> import sasmol.topology as topology
        >>> test = topology.CharmmTopology()
        >>> l1 = [1, 2, 4]
        >>> l2 = [2, 1, 4]
        >>> test.compare_list_ignore_order(l1, l2)
        True 

        Returns
        -------
        Boolean

        '''

        if len(l1) != len(l2):  # make sure the lenght of two lists are the same
            return False
        for item in l1:
            if l2.count(item) != 1:  # Dont allow duplicate elements during compare
                return False
        return True

    def check_charmm_atomic_order_reorganize(self, **kwargs):
        '''
        re-organize the atomic list according to the charmm topology contract
        do nothing if the atomic order alreay match that in the charmm topology file
        meanwhile make sure there are no missing or extra atoms
        H-atoms are required
        Patch N-ter for the first residue and C-ter for the last residue in each segment

        Parameters
        ----------
        
        kwargs 
            optional future arguments

        Examples
        -------
        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.check_charmm_atomic_order_reorganize()

        Returns
        -------
        error
            list : with error code

        Note
        ----
    
        error code needs to be addressed and handled globally

        '''
        
        error = []
        bin_path = config.__bin_path__
        self.read_charmm_topology(topology_file_path=bin_path + '/toppar/')
        self.setup_charmm_residue_atoms()
        self.initialize_children()
        children_segname = self.init_child('segnames')
        for child_segname in children_segname:
            child_segname.initialize_children()
            children = child_segname.init_child('resids')
            resid_nter = child_segname.resid()[0]
            resid_cter = child_segname.resid()[-1]
            for child in children:
                child_resname = child.resname()[0]
                child_resid = child.resid()[0]
                child_names = child.name()
                child_indices = child.index()
                #
                # print child_resname,self.charmm_residue_atoms[child_resname], child_names
                # print
                # self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname],
                # child_names)
                if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                    if child_resname == 'CYS':
                        child_resname = 'DISU'  # OPEN: try DISU if CYS doesnt work
                    if child_resname == 'HIS':
                        child_resname = 'HSE'  # OPEN: try HSE if HIS doesnt work
                        if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                            child_resname = 'HSD'  # OPEN: try HSD if HIS doesnt work
                            if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                                child_resname = 'HSP'  # OPEN: try HSP if HIS doesnt work
                    if child_resid == resid_nter:
                            # print child_resname
                        if child_resname == 'GLY':
                            patch = 'GLYP'
                        elif child_resname == 'PRO':
                            patch = 'PROP'
                        elif child_resname == 'PROP':
                            patch = 'PROP'
                        else:
                            patch = 'NTER'
                        self.patch_charmm_residue_atoms(child_resname, patch)
                        child_resname = child_resname + '_' + patch
                    if child_resid == resid_cter:
                        self.patch_charmm_residue_atoms(child_resname, 'CTER')
                        child_resname = child_resname + '_CTER'
                    if not self.compare_list_ignore_order(self.charmm_residue_atoms[child_resname], child_names):
                        error.append("For residue: " + child_resname + "\nthe atom names doesn't match those in the charmm topology file!\n" + str(
                            self.charmm_residue_atoms[child_resname]) + "\n" + str(child_names))
                        return error
                #
                # if self.charmm_residue_atoms[child_resname] == child_names:
                #	continue
                #import pprint
                # pprint.pprint(child_resname)
                # pprint.pprint(self.charmm_residue_atoms[child_resname],width=100)
                new_indices = []
                for name in self.charmm_residue_atoms[child_resname]:
                    new_indices.append(child_names.index(name))
                if len(child_indices) != len(new_indices):
                    error.append('Number of atoms doesnt match that in charmm topology file for \nresname: ' +
                                 child_resname + '\nand resid: ' + child_resid + '!')
                    return error
                for i in range(len(child_indices)):
                    # if child_indices[i]!=new_indices[i]:
                    self.atom()[child_indices[i] -
                                1] = child.atom()[new_indices[i]]
                    #self.index()[indices[i]-1] = child.index[new_indices[i]]
                    self.name()[child_indices[i] -
                                1] = child.name()[new_indices[i]]
                    self.loc()[child_indices[i] -
                               1] = child.loc()[new_indices[i]]
                    self.resname()[child_indices[i] -
                                   1] = child.resname()[new_indices[i]]
                    self.chain()[child_indices[i] -
                                 1] = child.chain()[new_indices[i]]
                    self.resid()[child_indices[i] -
                                 1] = child.resid()[new_indices[i]]
                    self.rescode()[child_indices[i] -
                                   1] = child.rescode()[new_indices[i]]
                    self.occupancy()[child_indices[i] -
                                     1] = child.occupancy()[new_indices[i]]
                    self.beta()[child_indices[i] -
                                1] = child.beta()[new_indices[i]]
                    self.segname()[child_indices[i] -
                                   1] = child.segname()[new_indices[i]]
                    self.element()[child_indices[i] -
                                   1] = child.element()[new_indices[i]]
                    self.charge()[child_indices[i] -
                                  1] = child.charge()[new_indices[i]]
                    # only frame-0 was handled in subset.init_child
                    self.coor()[0][child_indices[i] -
                                   1] = child.coor()[0][new_indices[i]]
        return error


