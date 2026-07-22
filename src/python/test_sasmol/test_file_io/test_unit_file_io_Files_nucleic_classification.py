'''Coordinate-only regression tests for the authoritative nucleic classifier.'''

import unittest

import numpy

import sasmol.system as system


class TestNucleicClassification(unittest.TestCase):

    def make_molecule(self, residues):
        '''Build residues as (resname, atom_names, resid, segname, chain, record).'''
        atom = []
        name = []
        resname = []
        resid = []
        segname = []
        chain = []
        record = []
        for this_resname, atom_names, this_resid, this_segname, this_chain, this_record in residues:
            for atom_name in atom_names:
                atom.append(this_record)
                name.append(atom_name)
                resname.append(this_resname)
                resid.append(this_resid)
                segname.append(this_segname)
                chain.append(this_chain)
                record.append(this_record)
        return system.Molecule_Maker(
            len(name), atom=record, name=name, resname=resname, resid=resid,
            segname=segname, chain=chain, moltype=['other'] * len(name))

    def one_residue_report(self, resname, atoms, record='ATOM'):
        molecule = self.make_molecule([(resname, atoms, 1, 'S', 'A', record)])
        report = molecule.classify_nucleic_moltypes()
        return molecule, report['segments']['S']

    def test_registry_is_mutually_exclusive_and_complete(self):
        molecule = system.Molecule()
        registry = molecule.nucleic_evidence_registry()
        classes = [registry['explicit_deoxy'], registry['explicit_ribonucleotide'],
                   registry['sugar_dependent_ambiguous'], registry['recognized_unspecified']]
        for index, one_class in enumerate(classes):
            for other_class in classes[index + 1:]:
                self.assertFalse(one_class.intersection(other_class))
        self.assertEqual(
            registry['recognized_nucleic'], set().union(*classes))
        for name in registry['recognized_nucleic']:
            self.assertNotEqual(molecule.moltype_for_resname(name), 'other')
        self.assertEqual(molecule.moltype_for_resname('NOT_A_RESIDUE'), 'other')
        self.assertEqual(molecule.moltype_for_resname('ALA'), 'protein')
        self.assertEqual(molecule.moltype_for_resname('TIP3'), 'water')

    def test_atom_aliases_are_exact(self):
        molecule = system.Molecule()
        self.assertEqual(molecule._normalize_nucleic_atom_name(" O2' "), "O2'")
        self.assertEqual(molecule._normalize_nucleic_atom_name(' O2* '), 'O2*')
        aliases = set(molecule.nucleic_rna_atom_aliases())
        for atom_name in ('O2', 'O2P', "O2''", 'O2**', "XO2'", "O2'X", ''):
            self.assertNotIn(molecule._normalize_nucleic_atom_name(atom_name), aliases)

    def test_residue_local_decision_table(self):
        complete_prime = ["C1'", "C2'", "C3'", "C4'", "O4'"]
        cases = [
            ('A', ["O2'"], 'rna', 'o2prime'),
            ('A', ['O2*'], 'rna', 'o2prime'),
            ('ADE', ["O2'"], 'rna', 'o2prime'),
            ('CYT', ['O2*'], 'rna', 'o2prime'),
            ('A', complete_prime, 'dna', 'complete_deoxy_sugar'),
            ('ADE', complete_prime, 'dna', 'complete_deoxy_sugar'),
            ('DA', ['P'], 'dna', 'explicit_deoxy_name'),
            ('DT', ['N3'], 'dna', 'explicit_deoxy_name'),
            ('A', ['N9'], 'ambiguous', 'insufficient_coordinate_evidence'),
            ('URA', ['N1', 'O2'], 'ambiguous', 'insufficient_coordinate_evidence'),
            ('THY', ['C5'], 'ambiguous', 'insufficient_coordinate_evidence'),
            ('DA', ["O2'"], 'conflict', 'conflict'),
            ('RNUA', complete_prime, 'conflict', 'conflict'),
        ]
        for resname, atoms, expected_identity, expected_evidence in cases:
            molecule, report = self.one_residue_report(resname, atoms)
            decision = report['residue_decisions'][0]
            self.assertEqual(decision['identity'], expected_identity, resname)
            self.assertEqual(decision['evidence'], expected_evidence, resname)
            self.assertEqual(list(molecule.moltype()),
                             ['nucleic' if expected_identity in ('ambiguous', 'conflict')
                              else expected_identity] * molecule.natoms())

    def test_incomplete_sugar_never_implies_dna(self):
        complete = ["C1'", "C2'", "C3'", "C4'", "O4'"]
        for missing_atom in complete:
            atoms = [atom for atom in complete if atom != missing_atom]
            molecule, report = self.one_residue_report('A', atoms)
            self.assertEqual(report['identity'], 'ambiguous')
            self.assertEqual(list(molecule.moltype()), ['nucleic'] * molecule.natoms())
        molecule, report = self.one_residue_report('A', ['P', 'O1P', 'O2P', "O3'", "O5'", "C5'"])
        self.assertEqual(report['identity'], 'ambiguous')

    def test_group_aggregation_canonicalizes_ambiguous_atoms(self):
        molecule = self.make_molecule([
            ('DA', ['P'], 1, 'D', 'A', 'ATOM'),
            ('ADE', ['N9'], 2, 'D', 'A', 'ATOM'),
            ('GUA', ['N7'], 3, 'D', 'A', 'HETATM')])
        report = molecule.classify_nucleic_moltypes()
        self.assertEqual(report['segments']['D']['identity'], 'resolved_dna')
        self.assertEqual(list(molecule.moltype()), ['dna', 'dna', 'dna'])
        self.assertEqual(report['segments']['D']['record_classes'], ['ATOM', 'HETATM'])

    def test_conflict_canonicalizes_every_recognized_atom(self):
        molecule = self.make_molecule([
            ('DA', ['P'], 1, 'X', 'A', 'ATOM'),
            ('A', ["O2'"], 2, 'X', 'A', 'ATOM'),
            ('GUA', ['N7'], 3, 'X', 'A', 'ATOM')])
        report = molecule.classify_nucleic_moltypes()
        self.assertEqual(report['segments']['X']['identity'], 'conflict')
        self.assertEqual(list(molecule.moltype()), ['nucleic', 'nucleic', 'nucleic'])

    def test_blank_identifiers_do_not_propagate(self):
        molecule = self.make_molecule([
            ('A', ["O2'"], 1, '', '', 'ATOM'),
            ('GUA', ['N7'], 2, '', '', 'ATOM')])
        report = molecule.classify_nucleic_moltypes()
        self.assertEqual(list(molecule.moltype()), ['rna', 'nucleic'])
        self.assertEqual(sorted(report['segments']), ['unassigned:1', 'unassigned:2'])
        self.assertEqual(report['segments']['unassigned:1']['grouping_source'], 'none')

    def test_reclassification_is_idempotent_and_uses_current_segnames(self):
        molecule = self.make_molecule([
            ('DA', ['P'], 1, 'X', 'A', 'ATOM'),
            ('A', ["O2'"], 2, 'X', 'B', 'ATOM')])
        first = molecule.classify_nucleic_moltypes()
        self.assertEqual(first['segments']['X']['identity'], 'conflict')
        first_moltypes = list(molecule.moltype())
        second = molecule.classify_nucleic_moltypes()
        self.assertEqual(list(molecule.moltype()), first_moltypes)
        self.assertEqual(second, first)
        molecule.setSegname(['D', 'R'])
        third = molecule.classify_nucleic_moltypes()
        self.assertEqual(third['segments']['D']['identity'], 'resolved_dna')
        self.assertEqual(third['segments']['R']['identity'], 'resolved_rna')
        self.assertEqual(list(molecule.moltype()), ['dna', 'rna'])

    def test_classification_changes_only_moltype_descriptors(self):
        molecule = self.make_molecule([('A', ["O2'", "C1'"], 1, 'R', 'A', 'HETATM')])
        snapshot = {
            'atom': list(molecule.atom()), 'name': list(molecule.name()),
            'resname': list(molecule.resname()), 'resid': list(molecule.resid()),
            'chain': list(molecule.chain()), 'segname': list(molecule.segname()),
            'coor': numpy.array(molecule.coor(), copy=True)}
        molecule.classify_nucleic_moltypes()
        self.assertEqual(list(molecule.atom()), snapshot['atom'])
        self.assertEqual(list(molecule.name()), snapshot['name'])
        self.assertEqual(list(molecule.resname()), snapshot['resname'])
        self.assertEqual(list(molecule.resid()), snapshot['resid'])
        self.assertEqual(list(molecule.chain()), snapshot['chain'])
        self.assertEqual(list(molecule.segname()), snapshot['segname'])
        self.assertTrue(numpy.array_equal(molecule.coor(), snapshot['coor']))


if __name__ == '__main__':
    unittest.main()
