#from __future__ import absolute_import
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
#	UTIL
#
#	12/10/2009	--	initial coding				:	jc
#	11/24/2011	-- 	moved to seperate file      :	jc
#	12/27/2015	-- 	refactored for release      :   jc
#	08/19/2016	-- 	added doc strings           :   jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Util holds general methods for file naming, os differences,
        chemical formula parsing, etc.
'''
#import system as system
import string
import copy
import numpy

def find_unique(this_list):
    return list(numpy.unique(this_list))

class Copy_Using_Mask():

    @classmethod
    def from_sasmol(class_instance,  mask, **kwargs):
        ''' 

        Parameters
        __________
        class_instance
            system object

        mask
            integer list

        kwargs
            optional future keyword arguments

        Returns
        _______
            molecule
                new system object

        Examples
        --------
            >>> import sasmol.system as system 
            >>> import sasmol.util as utilties 
            >>> molecule = system.Molecule('hiv1_gag.pdb')
            >>> molecule.name()[0] 
            'N'
            >>> molecule.name()[4] 
            'CA'
            >>> molecule.name()[14] 
            'HB1'
            >>> mask = [0, 4, 14]
            >>> new_molecule = utilities.Copy_Using_Mask.from_sasmol(molecule, mask) 
            >>> new_molecule.name()
            ['N', 'CA', 'HB1']
            >>> new_molecule.mass()
            array([ 14.00672 ,  12.01078 ,   1.007947])

        Note
        ____

            if more attributes are added in system.Atom() then the key lists below need
                to be updated

            currently only list_keys and numpy_keys are returned

            int_keys would have to be recalculated based on mask

            short_keys would have to be re-initialized (init_children)

            may want to consider passing a list of specific attributes to extract if memory is an issue

            copies all frames (untested)

            why can't this method be in subset?

        ''' 
       
        list_keys = ['_residue_flag', '_occupancy', '_charge', '_atom', '_chain', '_segname', '_beta', '_loc', '_element', '_name', '_rescode', '_moltype', '_resname']

        numpy_keys = ['_original_index', '_original_resid', '_index', '_resid', '_mass', '_coor'] 
     
        short_keys = ['_resnames', '_resids', '_elements', '_segnames', '_betas', '_names', '_moltypes', '_occupancies' ]      
        
        int_keys = ['_number_of_chains', '_number_of_betas', '_number_of_resids', '_number_of_names', '_number_of_moltypes', '_number_of_resnames', '_number_of_segnames', '_number_of_elements', '_id', '_number_of_occupancies' ]  
       
        other_keys = ['_header', '_conect', '_debug']  
            
        new_dict = {}

        number_of_frames = len(class_instance.coor())
        mask_length = len(mask) 
        
        natoms = class_instance._natoms
        all_data = [[] for x in xrange(natoms)]
        for i in xrange(natoms):
            if i in mask:
                count = 0
                for key, value in class_instance.__dict__.iteritems():
                    if key in list_keys:
                        all_data[count].append(value[i])
                        count += 1
        count = 0

        for key, value in class_instance.__dict__.iteritems():
            if key in list_keys:
                new_dict[key] = all_data[count]
                count += 1

        molecule = system.Molecule()
        molecule.__dict__ = new_dict  
      
        molecule.setId(0)

        molecule.setNatoms(mask_length) 
        
        molecule.setOriginal_index(numpy.take(class_instance.original_index(),mask))        
        molecule.setOriginal_resid(numpy.take(class_instance.original_resid(),mask))        
        molecule.setIndex(numpy.take(class_instance.index(),mask))        
        molecule.setResid(numpy.take(class_instance.resid(),mask))
                
        molecule.setMass(numpy.take(class_instance.mass(),mask))        
       
        molecule.setCoor(numpy.take(class_instance.coor(),mask,axis=1))
            
        return molecule

def duplicate_molecule(molecule, number_of_duplicates):
    return [copy.deepcopy(molecule) for x in xrange(number_of_duplicates)]


NAME, NUM, LPAREN, RPAREN, EOS = range(5)
import re
_lexer = re.compile(r"[A-Z][a-z]*|\d+|[()]|<EOS>").match
del re

# symbol, name, atomic number, molecular weight
_data = r"""'Ac', 'Actinium', 89, 227
'Ag', 'Silver', 47, 107.868
'Al', 'Aluminum', 13, 26.98154
'Am', 'Americium', 95, 243
'Ar', 'Argon', 18, 39.948
'As', 'Arsenic', 33, 74.9216
'At', 'Astatine', 85, 210
'Au', 'Gold', 79, 196.9665
'B', 'Boron', 5, 10.81
'Ba', 'Barium', 56, 137.33
'Be', 'Beryllium', 4, 9.01218
'Bi', 'Bismuth', 83, 208.9804
'Bk', 'Berkelium', 97, 247
'Br', 'Bromine', 35, 79.904
'C', 'Carbon', 6, 12.011
'Ca', 'Calcium', 20, 40.08
'Cd', 'Cadmium', 48, 112.41
'Ce', 'Cerium', 58, 140.12
'Cf', 'Californium', 98, 251
'Cl', 'Chlorine', 17, 35.453
'Cm', 'Curium', 96, 247
'Co', 'Cobalt', 27, 58.9332
'Cr', 'Chromium', 24, 51.996
'Cs', 'Cesium', 55, 132.9054
'Cu', 'Copper', 29, 63.546
'Dy', 'Dysprosium', 66, 162.50
'Er', 'Erbium', 68, 167.26
'Es', 'Einsteinium', 99, 254
'Eu', 'Europium', 63, 151.96
'F', 'Fluorine', 9, 18.998403
'Fe', 'Iron', 26, 55.847
'Fm', 'Fermium', 100, 257
'Fr', 'Francium', 87, 223
'Ga', 'Gallium', 31, 69.735
'Gd', 'Gadolinium', 64, 157.25
'Ge', 'Germanium', 32, 72.59
'H', 'Hydrogen', 1, 1.0079
'He', 'Helium', 2, 4.0026
'Hf', 'Hafnium', 72, 178.49
'Hg', 'Mercury', 80, 200.59
'Ho', 'Holmium', 67, 164.9304
'I', 'Iodine', 53, 126.9045
'In', 'Indium', 49, 114.82
'Ir', 'Iridium', 77, 192.22
'K', 'Potassium', 19, 39.0983
'Kr', 'Krypton', 36, 83.80
'La', 'Lanthanum', 57, 138.9055
'Li', 'Lithium', 3, 6.94
'Lr', 'Lawrencium', 103, 260
'Lu', 'Lutetium', 71, 174.96
'Md', 'Mendelevium', 101, 258
'Mg', 'Magnesium', 12, 24.305
'Mn', 'Manganese', 25, 54.9380
'Mo', 'Molybdenum', 42, 95.94
'N', 'Nitrogen', 7, 14.0067
'Na', 'Sodium', 11, 22.98977
'Nb', 'Niobium', 41, 92.9064
'Nd', 'Neodymium', 60, 144.24
'Ne', 'Neon', 10, 20.17
'Ni', 'Nickel', 28, 58.71
'No', 'Nobelium', 102, 259
'Np', 'Neptunium', 93, 237.0482
'O', 'Oxygen', 8, 15.9994
'Os', 'Osmium', 76, 190.2
'P', 'Phosphorous', 15, 30.97376
'Pa', 'Proactinium', 91, 231.0359
'Pb', 'Lead', 82, 207.2
'Pd', 'Palladium', 46, 106.4
'Pm', 'Promethium', 61, 145
'Po', 'Polonium', 84, 209
'Pr', 'Praseodymium', 59, 140.9077
'Pt', 'Platinum', 78, 195.09
'Pu', 'Plutonium', 94, 244
'Ra', 'Radium', 88, 226.0254
'Rb', 'Rubidium', 37, 85.467
'Re', 'Rhenium', 75, 186.207
'Rh', 'Rhodium', 45, 102.9055
'Rn', 'Radon', 86, 222
'Ru', 'Ruthenium', 44, 101.07
'S', 'Sulfur', 16, 32.06
'Sb', 'Antimony', 51, 121.75
'Sc', 'Scandium', 21, 44.9559
'Se', 'Selenium', 34, 78.96
'Si', 'Silicon', 14, 28.0855
'Sm', 'Samarium', 62, 150.4
'Sn', 'Tin', 50, 118.69
'Sr', 'Strontium', 38, 87.62
'Ta', 'Tantalum', 73, 180.947
'Tb', 'Terbium', 65, 158.9254
'Tc', 'Technetium', 43, 98.9062
'Te', 'Tellurium', 52, 127.60
'Th', 'Thorium', 90, 232.0381
'Ti', 'Titanium', 22, 47.90
'Tl', 'Thallium', 81, 204.37
'Tm', 'Thulium', 69, 168.9342
'U', 'Uranium', 92, 238.029
'Unh', 'Unnilhexium', 106, 263
'Unp', 'Unnilpentium', 105, 260
'Unq', 'Unnilquadium', 104, 260
'Uns', 'Unnilseptium', 107, 262
'V', 'Vanadium', 23, 50.9415
'W', 'Tungsten', 74, 183.85
'Xe', 'Xenon', 54, 131.30
'Y', 'Yttrium', 39, 88.9059
'Yb', 'Ytterbium', 70, 173.04
'Zn', 'Zinc', 30, 65.38
'Zr', 'Zirconium', 40, 91.22
'D', 'Deuterium', 1, 2.158"""


class Element:

    def __init__(self, symbol, name, atomicnumber, molweight):
        self.sym = symbol
        self.name = name
        self.ano = atomicnumber
        self.mw = molweight

    def getweight(self):
        return self.mw

    def addsyms(self, weight, result):
        result[self.sym] = result.get(self.sym, 0) + weight


def build_dict(s):
    import string
    answer = {}
    for line in string.split(s, "\n"):
        symbol, name, num, weight = eval(line)
        answer[symbol] = Element(symbol, name, num, weight)
    return answer


class ElementSequence:

    def __init__(self, *seq):
        self.seq = list(seq)
        self.count = 1

    def append(self, thing):
        self.seq.append(thing)

    def getweight(self):
        sum = 0.0
        for thing in self.seq:
            sum = sum + thing.getweight()
        return sum * self.count

    def set_count(self, n):
        self.count = n

    def __len__(self):
        return len(self.seq)

    def addsyms(self, weight, result):
        totalweight = weight * self.count
        for thing in self.seq:
            thing.addsyms(totalweight, result)

    def displaysyms(self, sym2elt):
        result = {}
        self.addsyms(1, result)
        items = result.items()
        items.sort()
        for sym, count in items:
            print(sym, " : ", count)


class Tokenizer:

    def __init__(self, input):
        self.input = input + "<EOS>"
        self.i = 0

    def gettoken(self):
        global ttype, tvalue
        self.lasti = self.i
        m = _lexer(self.input, self.i)
        if m is None:
            self.error("unexpected character")
        self.i = m.end()
        tvalue = m.group()
        if tvalue == "(":
            ttype = LPAREN
        elif tvalue == ")":
            ttype = RPAREN
        elif tvalue == "<EOS>":
            ttype = EOS
        elif "0" <= tvalue[0] <= "9":
            ttype = NUM
            tvalue = int(tvalue)
        else:
            ttype = NAME

    def error(self, msg):
        emsg = msg + ":\n"
        emsg = emsg + self.input[:-5] + "\n"  # strip <EOS>
        emsg = emsg + " " * self.lasti + "^\n"
        raise ValueError(emsg)


def parse(s, sym2elt):
    global t, ttype, tvalue
    t = Tokenizer(s)
    t.gettoken()
    seq = parse_sequence(sym2elt)
    if ttype != EOS:
        t.error("expected end of input")
    return seq


def parse_sequence(sym2elt):
    global t, ttype, tvalue
    seq = ElementSequence()
    while ttype in (LPAREN, NAME):
        # parenthesized expression or straight name
        if ttype == LPAREN:
            t.gettoken()
            thisguy = parse_sequence(sym2elt)
            if ttype != RPAREN:
                t.error("expected right paren")
            t.gettoken()
        else:
            assert ttype == NAME
            if sym2elt.has_key(tvalue):
                thisguy = ElementSequence(sym2elt[tvalue])
            else:
                t.error("'" + tvalue + "' is not an element symbol")
            t.gettoken()
    # followed by optional count
        if ttype == NUM:
            thisguy.set_count(tvalue)
            t.gettoken()
        seq.append(thisguy)
    if len(seq) == 0:
        t.error("empty sequence")
    return seq

def get_chemical_formula(formula_string):

    Atomic = properties.Atomic()
    #standard_atomic_weights = Atomic.amu(keep_lower_case=True)
    amu = Atomic.amu(keep_lower_case=True)
    sym2elt = build_dict(_data)

#    #sym2elt = amu

    formula_dictionary = {}

    error = []
    try:
        seq = parse(formula_string.strip(" "),sym2elt)
        #seq.displaysyms(sym2elt)
        seq.addsyms(1,formula_dictionary)
        items = formula_dictionary.items()
        items.sort()

        #for sym, count in items:
        #    print sym," :: ",count

    except ValueError, detail:
        print(str(detail))
        error.append(detail)

    return error,formula_dictionary

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

if __name__ == "__main__":

    path = './'

    filename = 'run_0'

    '''
    while 1:
        x = raw_input("? ")
        fields = string.split(x)
        if len(fields) != 2:
            print("input must have two fields")
            continue
        action, formula = fields
        ok = 0
        try:
            seq = parse(formula,sym2elt)
            ok = 1
        except ValueError, detail:
            print(str(detail))
        if not ok:
            continue
        if action == "molw":
            print("molecular weight", seq.getweight())
        elif action == "syms":
            seq.displaysyms(sym2elt)
        else:
            print("unknown action:", action)
    '''
