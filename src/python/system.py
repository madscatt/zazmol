# from __future__ import unicode_literals
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
#   SYSTEM
#
# 12/4/2009	--	initial coding			:	jc
# 12/10/2009	--	doc strings 			:	jc
# 01/11/2010	--	new design pattern		:	jc
# 12/25/2015	--	refactored for release  :   jc
# 07/23/2016	--	refactored for Python 3 :   jc
# 08/19/2016	--	added doc strings       :   jc
#
# 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
# *      **
'''
	System is the main module that contains the base classes that
	describe the objects in the system.  An instance of the system
	can be created of Atom, Molecule, or System classes.

	For purely atomic based systems, a set of utilities are provided
	to hold information, to calculate properties, and to manipulate
	the structures in space.  In addition, basic file input / output
	is provided in the Mol class (and inherited up the classes
	as dictated by the calling code).

'''

import sys
import os
import numpy
import sasmol.file_io as file_io
import sasmol.calculate as calculate
import sasmol.operate as operate
import sasmol.subset as subset
import sasmol.properties as properties
import sasmol.topology as topology
import sasmol.view as view

import sasmol.config as config


class Atom(
        file_io.Files,
        calculate.Calculate,
        operate.Move,
        subset.Mask,
        properties.Atomic,
        topology.Topology,
        view.View):

    """Base class for molecular system objects.

    ``Atom`` mixes in file I/O, coordinate operations, subset/mask operations,
    properties, topology, and viewing helpers used across SASMOL.

    Accessor pattern
    ----------------
    Most fields are exposed in getter/setter pairs:
    ``field_name()`` returns the current value and ``setField_name(value)``
    (or legacy camel-case variants such as ``setName``) assigns a new value.
    Many of these values are atom-aligned arrays/lists with length ``natoms``.

    Parameters
    ----------
    args
        optional integer : object id

    kwargs
        optional keyword arguments:
        ``filename`` (string), ``id`` (integer), ``debug`` (boolean)

    Returns
    -------
    system object
        if called with ``filename`` (kwarg or positional string), returns
        an initialized system object with data read using ``read_pdb()``.

    Examples
    --------
    >>> import sasmol.system as system
    >>> molecule = system.Molecule(filename='hiv1_gag.pdb')
    >>> molecule.index()[0]
    1
    >>> molecule.setIndex([x - 10 for x in molecule.index()])
    >>> molecule.index()[0]
    -9

    Notes
    -----
    ``self`` is omitted from the parameter section in generated docs.

    """

    def __init__(self, *args, **kwargs):
        self._filename = kwargs.pop('filename', None)
        self._id = kwargs.pop('id', 0)
        self._debug = kwargs.pop('debug', None)

        self._total_mass = 0.0
        self._natoms = 0
        self._mass = None
        self._coor = None
        self._center_of_mass = None
        self._conect = []

        if config.__logging_level__ == 'DEBUG':
            self._debug = True

        self._defined_with_input_file = False
        self._argument_flag = False
        self._id_flag = False

        try:
            if self._filename is not None:
                if (os.path.isfile(self._filename)):
                    self.read_pdb(self._filename)
                    self._defined_with_input_file = True
            else:

                for argument in args:
                    self._argument_flag = True
                    if (os.path.isfile(str(argument))):
                        try:
                            self.read_pdb(argument)
                            self._defined_with_input_file = True
                            self._filename = argument
                            break
                        except Exception:
                            pass
                    else:
                        try:
                            self._id = int(argument)
                            self._id_flag = True
                            break
                        except Exception:
                            pass

        except Exception:
            pass

    def __repr__(self):

        if self._defined_with_input_file:
            return "sasmol object initialied with filename = " + self._filename
        elif self._argument_flag and not self._id_flag:
            return "sasmol object (no file found)"
        else:
            return "sasmol object"

    def setFilename(self, newValue):
        self._filename = newValue

    def filename(self):
        return self._filename

    def __add__(self, other, **kwargs):
        '''
        Combine another molecule into this molecule in-place.

        Parameters
        ----------
        other
            system object : molecule-like object to append

        Returns
        -------
        None
            the current object is modified in-place

        Examples
        --------
        >>> import sasmol.system as system
        >>> molecule_1 = system.Molecule(filename='hiv1_gag.pdb')
        >>> molecule_2 = system.Molecule(filename='lysozyme.pdb')
        >>> molecule_1 + molecule_2
        >>> molecule_1.natoms() > 0
        True

        Notes
        -----
        Descriptors missing from ``other`` are skipped.
        ``_natoms`` is recomputed from ``_name`` and ``_index`` is re-sequenced.
        Frame-count compatibility is not validated here.
        '''

        for key, value in self.__dict__.items():
            try:
                if type(value) is list:
                    self.__dict__[key].extend(other.__dict__[key])
                elif type(value) is numpy.ndarray:
                   # other.__dict__[key]
                    if key == '_coor':
                        self.__dict__[key] = numpy.concatenate(
                            (self.__dict__[key], other.__dict__[key]), axis=1)
                    elif len(value.shape) == 1:
                        self.__dict__[key] = numpy.concatenate(
                            (self.__dict__[key], other.__dict__[key]))
                    else:
                        print(
                            'numpy array not added for self.__dict__[key]: ' +
                            str(key))

            except Exception:
                pass

            self._natoms = len(self._name)
            self._index = numpy.array(
                [x + 1 for x in range(self._natoms)], numpy.int32)

    def setId(self, newValue):
        self._id = newValue

    def id(self):
        return self._id

        # properties

    def number_of_frames(self):
        return self._coor[:, 0, 0].shape[0]

    def setNumber_of_frames(self, newValue):
        self._number_of_frames = newValue

    def atom(self):
        return self._atom

    def setAtom(self, newValue):
        self._atom = newValue

    def index(self):
        return self._index

    def setIndex(self, newValue):
        self._index = newValue

    def original_index(self):
        return self._original_index

    def setOriginal_index(self, newValue):
        self._original_index = newValue

    def original_resid(self):
        return self._original_resid

    def setOriginal_resid(self, newValue):
        self._original_resid = newValue

    def name(self):
        return self._name

    def setName(self, newValue):
        self._name = newValue

    def loc(self):
        return self._loc

    def setLoc(self, newValue):
        self._loc = newValue

    def resname(self):
        return self._resname

    def setResname(self, newValue):
        self._resname = newValue

    def chain(self):
        return self._chain

    def setChain(self, newValue):
        self._chain = newValue

    def resid(self):
        return self._resid

    def setResid(self, newValue):
        self._resid = newValue

    def coor(self):
        return self._coor

    def setCoor(self, newValue):
        self._coor = newValue

    def rescode(self):
        return self._rescode

    def setRescode(self, newValue):
        self._rescode = newValue

    def occupancy(self):
        return self._occupancy

    def setOccupancy(self, newValue):
        self._occupancy = newValue

    def beta(self):
        return self._beta

    def setBeta(self, newValue):
        self._beta = newValue

    def segname(self):
        return self._segname

    def setSegname(self, newValue):
        self._segname = newValue

    def element(self):
        return self._element

    def setElement(self, newValue):
        self._element = newValue

    def charge(self):
        return self._charge

    def setCharge(self, newValue):
        self._charge = newValue

    def atom_charge(self):
        return self._atom_charge

    def setAtom_charge(self, newValue):
        self._atom_charge = newValue

    def atom_vdw(self):
        return self._atom_vdw

    def setAtom_vdw(self, newValue):
        self._atom_vdw = newValue

    def formula(self):
        return self._formula

    def setFormula(self, newValue):
        self._formula = newValue

    def mass(self):
        return self._mass

    def setMass(self, newValue):
        self._mass = newValue

    def total_mass(self):
        if (self._total_mass is None):
            self._total_mass = calculate.Calculate.calculate_mass(self)
        return self._total_mass

    def setTotal_mass(self, newValue):
        self._total_mass = newValue

    def natoms(self):
        return self._natoms

    def setNatoms(self, newValue):
        self._natoms = newValue

    def number_of_atoms(self):
        return self._number_of_atoms

    def setNumber_of_atoms(self, newValue):
        self._number_of_atoms = newValue

    def moltype(self):
        return self._moltype

    def setMoltype(self, newValue):
        self._moltype = newValue

    # --- diagnostics ---
    def check_integrity(self, fast_check=False):
        import sasmol.utilities as utilities
        return utilities.check_integrity(self, fast_check=fast_check)

    def validate_integrity(self, fast_check=False):
        '''
        Validate that core atom-aligned descriptors match ``natoms``.

        Parameters
        ----------
        fast_check
            bool : reserved flag forwarded to integrity checks

        Returns
        -------
        dict
            mapping of descriptor names to observed lengths

        Raises
        ------
        ValueError
            if required descriptors are missing or have mismatched lengths

        Examples
        --------
        >>> import sasmol.system as system
        >>> molecule = system.Molecule_Maker(2, name=['N', 'CA'])
        >>> result = molecule.validate_integrity()
        >>> 'name' in result
        True
        '''
        import sasmol.utilities as utilities

        results = utilities.check_integrity(
            self, fast_check=fast_check, warn=False)
        natoms = self.natoms()
        errors = []

        for key, length in results.items():
            if length is None:
                errors.append('%s is missing' % key)
            elif length == 'N/A':
                errors.append('%s has no length' % key)
            elif length != natoms:
                errors.append('%s has length %s, expected %s' %
                              (key, length, natoms))

        if errors:
            raise ValueError('Integrity check failed: ' + '; '.join(errors))

        return results

    def record(self):
        self._moltype = newValue

    def setRecord(self, newValue):
        self._atom = newValue

    def altloc(self):
        return self._loc

    def setAltloc(self, newValue):
        self._loc = newValue

    def icode(self):
        return self._rescode

    def setIcode(self, newValue):
        self._rescode = newValue


class Molecule(Atom):

    """Molecule class describing one molecular object.

    ``Molecule`` inherits all functionality from ``Atom`` and is the primary
    entry point used by callers.

    Parameters
    ----------
    args
        optional integer : object id

    kwargs
        optional keyword arguments:
        ``filename`` (string), ``id`` (integer), ``debug`` (boolean)

    Returns
    -------
    system object
        molecule instance, optionally initialized from a PDB file

    Examples
    --------
    >>> import sasmol.system as system
    >>> molecule = system.Molecule(filename='hiv1_gag.pdb')
    >>> empty = system.Molecule()
    >>> identified = system.Molecule(id=7)
    """

    def __init__(self, *args, **kwargs):
        Atom.__init__(self, *args, **kwargs)

    def fasta(self):
        return self._fasta

    def setFasta(self, newValue):
        self._fasta = newValue

    def setCharmm_type(self, newValue):
        self._charmm_type = newValue

    def charmm_type(self):
        return self._charmm_type

    def header(self):
        return self._header

    def setHeader(self, newValue):
        self._header = newValue

    def unitcell(self):
        return self._unitcell

    def setUnitcell(self, newValue):
        self._unitcell = newValue

    def rg(self):
        return self._rg

    def setRg(self, newValue):
        self._rg = newValue

    def pmi(self):
        return self._pmi

    def setPmi(self, newValue):
        self._pmi = newValue

    def minimum(self):
        return self._minimum

    def setMinimum(self, newValue):
        self._minimum = newValue

    def maximum(self):
        return self._maximum

    def setMaximum(self, newValue):
        self._maximum = newValue

    def one_letter_resname(self):
        return self._one_letter_resname

    def setOne_letter_resname(self, newValue):
        self._one_letter_resname = newValue

    def center_of_mass(self):
        return self._center_of_mass

    def setCenter_of_mass(self, newValue):
        self._center_of_mass = newValue

    def names(self):
        return self._names

    def setNames(self, newValue):
        self._names = newValue

    def resnames(self):
        return self._resnames

    def setResnames(self, newValue):
        self._resnames = newValue

    def resids(self):
        return self._resids

    def setResids(self, newValue):
        self._resids = newValue

    def chains(self):
        return self._chains

    def setChains(self, newValue):
        self._chains = newValue

    def segnames(self):
        return self._segnames

    def setSegnames(self, newValue):
        self._segnames = newValue

    def occupancies(self):
        return self._occupancies

    def setOccupancies(self, newValue):
        self._occupancies = newValue

    def betas(self):
        return self._betas

    def setBetas(self, newValue):
        self._betas = newValue

    def moltypes(self):
        return self._moltypes

    def setMoltypes(self, newValue):
        self._moltypes = newValue

    def elements(self):
        return self._elements

    def setElements(self, newValue):
        self._elements = newValue

    def names_mask(self):
        return self._names_mask

    def setNames_mask(self, newValue):
        self._names_mask = newValue

    def resnames_mask(self):
        return self._resnames_mask

    def setResnames_mask(self, newValue):
        self._resnames_mask = newValue

    def resids_mask(self):
        return self._resids_mask

    def setResids_mask(self, newValue):
        self._resids_mask = newValue

    def chains_mask(self):
        return self._chains_mask

    def setChains_mask(self, newValue):
        self._chains_mask = newValue

    def occupancies_mask(self):
        return self._occupancies_mask

    def setOccupancies_mask(self, newValue):
        self._occupancies_mask = newValue

    def betas_mask(self):
        return self._betas_mask

    def setBetas_mask(self, newValue):
        self._betas_mask = newValue

    def elements_mask(self):
        return self._elements_mask

    def setElements_mask(self, newValue):
        self._elements_mask = newValue

    def segnames_mask(self):
        return self._segnames_mask

    def setSegnames_mask(self, newValue):
        self._segnames_mask = newValue

    def number_of_names(self):
        return self._number_of_names

    def setNumber_of_names(self, newValue):
        self._number_of_names = newValue

    def number_of_resnames(self):
        return self._number_of_resnames

    def setNumber_of_resnames(self, newValue):
        self._number_of_resnames = newValue

    def number_of_resids(self):
        return self._number_of_resids

    def setNumber_of_resids(self, newValue):
        self._number_of_resids = newValue

    def number_of_chains(self):
        return self._number_of_chains

    def setNumber_of_chains(self, newValue):
        self._number_of_chains = newValue

    def number_of_segnames(self):
        return self._number_of_segnames

    def setNumber_of_segnames(self, newValue):
        self._number_of_segnames = newValue

    def number_of_occupancies(self):
        return self._number_of_occupancies

    def setNumber_of_occupancies(self, newValue):
        self._number_of_occupancies = newValue

    def number_of_betas(self):
        return self._number_of_betas

    def setNumber_of_betas(self, newValue):
        self._number_of_betas = newValue

    def number_of_moltypes(self):
        return self._number_of_moltypes

    def setNumber_of_moltypes(self, newValue):
        self._number_of_moltypes = newValue

    def number_of_elements(self):
        return self._number_of_elements

    def setNumber_of_elements(self, newValue):
        self._number_of_elements = newValue

    def residue_charge(self):
        return self._residue_charge

    def setResidue_charge(self, newValue):
        self._residue_charge = newValue

    def conect(self):
        return self._conect

    def setConect(self, newValue):
        self._conect = newValue

    def residue_flag(self):
        try:
            return self._residue_flag
        except Exception:
            self._residue_flag = False
            return self._residue_flag

    def setResidue_flag(self, newValue):
        self._residue_flag = newValue

    def set_average_vdw(self):
        import sasmol.operate as operate
        operate.set_average_vdw(self)
        return


class System(Atom):

    """ System is a class that is used to aggregate all components. It inherits
        all of attributes from Atom.

        Class has several initialization options

        Parameters
        ----------
        args
            optional integer : self._id

        kwargs
            optional keyword arguments
                string filename (filename = 'hiv1_gag.pdb') : default = None
                integet id (id=3) : default = 0
                boolean debug (debug = True) : default = None

        Returns
        -------
        system object
            if called with string (or with filename kwarg) returns
            an initialized system object with data read in using
            file_io.read_pdb()

        Examples
        -------

        >>> import sasmol.system as system
        >>> molecule = system.Molecule(filename='hiv1_gag.pdb')


        >>> molecule = system.Molecule()
        >>> molecule = system.Molecule(id=7)
        >>> molecule = system.Molecule(debug=True)
        >>> molecule = system.Molecule('hiv1_gag.pdb')
        >>> molecule = system.Molecule(filename='hiv1_gag.pdb', id=0, debug=False)

    """

    def __init__(self, *args, **kwargs):
        Atom.__init__(self, *args, **kwargs)


class Molecule_Maker(Atom):
    """
    Build a molecule from constructor-provided atom descriptors.

    ``Molecule_Maker`` is intended for programmatic molecule construction when
    no source PDB/DCD file is being read.

    For most atom fields, pass either:
    1) a scalar value applied to every atom, or
    2) a per-atom list/tuple/array with length exactly ``natoms``.

    Parameters
    ----------
    natoms
        integer : number of atoms in the molecule

    atom, name, loc, resname, chain, rescode, occupancy, beta, segname, element, charge, moltype
        scalar values or per-atom sequences

    index, resid
        integer sequence fields; default values are generated when omitted

    coor
        numpy float array : coordinate array

    kwargs
        optional future arguments

    Returns
    -------
    system object

    Examples
    --------
    >>> import sasmol.system as system
    >>> molecule = system.Molecule_Maker(4, name='Ar', segname='ARG0')
    >>> molecule = system.Molecule_Maker(2, name=['N', 'CA'], resid=[1, 1])
    >>> molecule.validate_integrity()

    """

    def _expand_constructor_value(self, field_name, value, natoms):
        if isinstance(value, (list, tuple, numpy.ndarray)):
            length = len(value)

            if length != natoms:
                raise ValueError(
                    '%s must be a scalar value or have length %d' %
                    (field_name, natoms))

            return list(value)

        return [value for x in range(natoms)]

    def _validate_constructor_length(self, field_name, value, natoms):
        if not isinstance(value, (list, tuple, numpy.ndarray)):
            raise ValueError('%s must have length %d' % (field_name, natoms))

        if len(value) != natoms:
            raise ValueError('%s must have length %d' % (field_name, natoms))

    def __init__(
            self,
            natoms,
            atom='ATOM',
            index=None,
            name='C',
            loc=' ',
            resname='DUM',
            chain='A',
            resid=None,
            rescode=' ',
            coor=None,
            occupancy='0.00',
            beta='0.00',
            segname='DUM',
            element='C',
            charge=' ',
            moltype='protein',
            residue_flag=False,
            **kwargs):

        Atom.__init__(self)

        self._natoms = natoms
        self._number_of_atoms = natoms

        self._atom = self._expand_constructor_value('atom', atom, natoms)

        if index is not None:
            self._validate_constructor_length('index', index, natoms)
            self._index = index
        else:
            self._index = numpy.array(
                [x + 1 for x in range(natoms)], numpy.int32)

        self._name = self._expand_constructor_value('name', name, natoms)
        self._loc = self._expand_constructor_value('loc', loc, natoms)
        self._resname = self._expand_constructor_value(
            'resname', resname, natoms)
        self._chain = self._expand_constructor_value('chain', chain, natoms)

        if resid is not None:
            self._validate_constructor_length('resid', resid, natoms)
            self._resid = resid
        else:
            self._resid = numpy.array(
                [x + 1 for x in range(natoms)], numpy.int32)

        self._rescode = self._expand_constructor_value(
            'rescode', rescode, natoms)

        if coor is not None:
            self._coor = coor
        else:
            self._coor = numpy.zeros((1, natoms, 3), config.COORD_DTYPE)

        self._occupancy = self._expand_constructor_value(
            'occupancy', occupancy, natoms)
        self._beta = self._expand_constructor_value('beta', beta, natoms)
        self._charge = self._expand_constructor_value(
            'charge', charge, natoms)
        self._segname = self._expand_constructor_value(
            'segname', segname, natoms)
        self._element = self._expand_constructor_value(
            'element', element, natoms)

        self.setOriginal_index(self.index())
        self.setOriginal_resid(self.resid())

        self._residue_flag = self._expand_constructor_value(
            'residue_flag', residue_flag, natoms)
        self._moltype = self._expand_constructor_value(
            'moltype', moltype, natoms)

        self._total_mass = calculate.Calculate.calculate_mass(self)

    def residue_flag(self):
        return self._residue_flag
