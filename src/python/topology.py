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

import numpy
import sasmol.charmm_topology as charmm_topology
import sasmol.properties as properties

class Topology(charmm_topology.CharmmTopology):

    '''
    This class contains topology information used other modules.

    '''

    def __init__(self):
        pass

    def renumber(self, **kwargs):
        '''
        Method to renumber index and resid fields

        Parameters
        ----------
        kwargs 
            index :
            resid :
                                                                                     
        Returns
        -------
        None
            updated system object 
                index
                resid

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.index()[0]
        1
        >>> molecule.resid()[0]
        1
       
        default renumber() with no arguments renumbers
        both index and resid starting at 1 
        
        >>> molecule.renumber()
        >>> molecule.index()[0]
        1
        >>> molecule.resid()[0]
        1

        passing a kwarg with a value will renumber appropriate
        attribute with value passed

        >>> molecule.renumber(index=3)
        >>> molecule.index()[0]
        3
        >>> molecule.resid()[0]
        1

        >>> molecule.renumber(resid=8)
        >>> molecule.index()[0]
        3
        >>> molecule.resid()[0]
        8

        >>> molecule.renumber(index=223, resid=18)
        >>> molecule.index()[0]
        223
        >>> molecule.resid()[0]
        18

        '''

        renumber_index_flag = False
        renumber_index_start = 1
       
        renumber_resid_flag = False
        renumber_resid_start = 1     
    
        try:
            if kwargs['index']:
                renumber_index_flag = True
                renumber_index_start = kwargs['index']
        except:
            pass

        try:
            if kwargs['resid']:
                renumber_resid_flag = True     
                renumber_resid_start = kwargs['resid']
        except:
            pass

        if renumber_index_flag:
            start = renumber_index_start
            index = numpy.array([x for x in xrange(start,start+self._natoms)])
            self._index = index

        if renumber_resid_flag:
            resid = self._resid
            resid_array=[] 
            count = renumber_resid_start
            for i in xrange(len(resid)):
                this_resid = resid[i]
                if(i==0):
                    last_resid = this_resid
                else:
                    if(this_resid != last_resid):	
                        count += 1
                resid_array.append(count)	
                last_resid = this_resid

            self._resid = numpy.array(resid_array, numpy.int)

        return

    def make_constraint_pdb(self, filename, basis_type, **kwargs):
        ''' 
     
        Method to rename attribute fields and assign values
        of 1.00 for the given basis type.

        "beta" is the default attribute field to reassign
        
        Default usage sets all attribute fields to zero
        prior to re-assigning basis type atoms to 1.00

        Parameters
        ----------
        filename    string : name of output PDB file to write 
             
        basis types: 

            'heavy : all atoms except hydrogen
          
            'protein' : all atoms in moltype protein
             
            'nucleic' : all atoms in moltype nucleic
                                        
            'backbone'
                -> proteins: N, CA, C, O
            
            'solute'    : all protein and nucleic


        kwargs 
            optional arguments
                
                field='beta' : write 1.00 to beta field

                field='occupancy' : write 1.00 to occupancy field

                reset=True : set all values to 0.00 before setting
                                 individual vaules to 1.00

                reset=False: only change desired values to 1.00
                                 ignoring all other values
                    
        Returns
        -------
        None
            updated system object
            writes pdb to disk

        Examples
        --------
        >>> import sasmol.system as system
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule.make_constraint_pdb('constrainted_backbone_hiv1_gag.pdb', 'backbone')

        >>> molecule.make_constraint_pdb('constrainted_backbone_occupancy_hiv1_gag.pdb', 'backbone', field='occupancy')

        constrtain all heavy atoms

        >>> molecule.make_constraint_pdb('constrainted_heavy_hiv1_gag.pdb', 'heavy')

        constrain all protein atoms
        
        >>> molecule.make_constraint_pdb('constrainted_protein_hiv1_gag.pdb', 'protein')

        constrain all non-solvent atoms
        
        >>> molecule.make_constraint_pdb('constrainted_solute_hiv1_gag.pdb', 'solute')
        
        Note
        ----

        Only moltype `protein` and `nucleic` are supported.

        Output PDB file will contain coordinates from frame 0
             
        '''

        try:
            field = kwargs['field']
        except:
            field = 'beta'

        try:
            reset = kwargs['reset']
        except:
            reset = True
    
        if reset:
            new_value = ['0.00' for x in xrange(self._natoms)]
            if field == 'beta':
                self._beta = new_value
            elif field == 'occupancy':
                self._occupancy = new_value

        if basis_type == 'backbone':
            basis_string = "name[i] == 'N' or name[i] == 'CA' or name[i] == 'C' or name[i] == 'O'"
                           
        elif basis_type == 'heavy':
            basis_string = "name[i][0] != 'H'"

        elif basis_type == 'protein':
            basis_string = "moltype[i] == 'protein'"
        
        elif basis_type == 'nucleic':
            basis_string = "moltype[i] == 'nucleic'"

        elif basis_type == 'solute':
            basis_string = "moltype[i] == 'protein' or moltype[i] == 'nucleic'"

        error, mask = self.get_subset_mask(basis_string)

        if len(error) > 0:
            print('error = ', error)
            return

        if field == 'beta':
            descriptor = self._beta
        elif field == 'occupancy':
            descriptor = self._occupancy

        error = self.set_descriptor_using_mask(mask, descriptor, '1.00') 
        
        if len(error) > 0:
            print('error = ', error)
            return
               
        if field == 'beta': 
            self._beta = descriptor
        elif field == 'occupancy': 
            self._occupancy = descriptor

        self.write_pdb(filename,0,'w')

        return


    def make_backbone_pdb_from_fasta(self, filename, moltype, **kwargs):
        
        """
        Method to write a PDB file of one atom per residue
        based on input FASTA sequence

        http://en.wikipedia.org/wiki/FASTA_format

        Parameters
        ----------
        filename    string : name of output PDB file to write 
        
        moltype string
                -> 'protein' ; fasta sequence of a protein 
                -> 'nucleic' ; fasta sequence of a nucleic acid 
             
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

        >>> molecule.make_backbone_pdb_from_fasta('sequence.pdb', 'protein')


        Note
        ____

        Only protein and nucleic acid one-letter codes are supported
        Protein PDB files are written with one "CA" atom per residue
        Nucleic acid PDF files are written with on "O5'" atom per residue
       
        N-terminal patches for proteins (GLYP, PROP) are accommodated
        
        All coordinate values are set to 0.000 
         
        """

        if moltype == "protein":
            residue_dictionary = self.one_to_three_letter_protein_residue_dictionary            
            sequence_name = "CA" 

        elif moltype == "nucleic":
            residue_dictionary = self.one_to_three_letter_nucleic_residue_dictionary 
            sequence_name = "O5'" 

        sequence = []
        first_flag = True

        for residue in self._fasta:

            if moltype == "protein":
                if residue == 'G' and first_flag:
                    sequence.append('GLYP')
                elif residue == 'P' and first_flag:
                    sequence.append('PROP')
                else:     
                    sequence.append(residue_dictionary[residue])

            elif moltype == "nucleic":
                sequence.append(residue_dictionary[residue])
        
            first_flag = False

        import sasmol.system as system
         
        molecule = system.Molecule_Maker(len(sequence), name=sequence_name)

        molecule._resname = sequence

        molecule.write_pdb(filename, 0, 'w')

        return

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
        >>> print(molecule.fasta()[:5]) # doctest: +NORMALIZE_WHITESPACE
        ['G', 'A', 'R', 'A', 'S'] 

        >>> molecule.create_fasta(fasta_format=True)
        >>> print(molecule.fasta()[:5])
        >
        GAR
        >>> molecule.create_fasta(fasta_format=True,width='60')
        >>> print(molecule.fasta()[:-1])
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
        >>> print(molecule.fasta()[:90])
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



if __name__ == "__main__":
    import doctest
    #doctest.testmod(verbose=True)
    doctest.testmod()




