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

import unittest

import sasmol.system as system


class Test_unit_file_io_Files_moltype_report(unittest.TestCase):

    def test_ambiguous_nucleic_resnames_are_reported_without_mutation(self):
        molecule = system.Molecule_Maker(
            3,
            name=['N9', 'N1', 'N7'],
            resname=['ADE', 'CYT', 'GUA'],
            resid=[1, 2, 3],
            segname='DNA1',
            moltype=['nucleic', 'nucleic', 'nucleic'])

        report = molecule.moltype_by_segname_report()

        self.assertEqual(report['overall_status'], 'ambiguous')
        self.assertEqual(molecule.moltype(), ['nucleic', 'nucleic', 'nucleic'])
        segment = report['segments']['DNA1']
        self.assertEqual(segment['status'], 'ambiguous')
        self.assertEqual(segment['assigned_moltypes'], ['nucleic'])
        self.assertEqual(segment['ambiguous_resnames'], ['ADE', 'CYT', 'GUA'])
        self.assertEqual(segment['residue_count'], 3)
        self.assertIn('insufficient coordinate evidence',
                      segment['evidence'][0])

    def test_specific_nucleic_resnames_keep_segment_clean(self):
        molecule = system.Molecule_Maker(
            2,
            name=['P', 'C5'],
            resname=['DA', 'DT'],
            resid=[1, 2],
            segname='DNA1',
            moltype=['dna', 'dna'])

        report = molecule.moltype_by_segname_report()

        self.assertEqual(report['overall_status'], 'clean')
        segment = report['segments']['DNA1']
        self.assertEqual(segment['status'], 'resolved_dna')
        self.assertEqual(segment['dna_resname_evidence'], ['DA', 'DT'])

    def test_rna_atom_name_evidence_can_resolve_overlap_resnames(self):
        molecule = system.Molecule_Maker(
            2,
            name=["O2'", 'N9'],
            resname=['ADE', 'GUA'],
            resid=[1, 2],
            segname='RNA1',
            moltype=['rna', 'rna'])

        report = molecule.moltype_by_segname_report()

        self.assertEqual(report['overall_status'], 'clean')
        segment = report['segments']['RNA1']
        self.assertEqual(segment['status'], 'resolved_rna')
        self.assertEqual(segment['ambiguous_resnames'], ['GUA'])
        self.assertEqual(segment['rna_atom_evidence'], ["O2'"])

    def test_mixed_segment_is_reported(self):
        molecule = system.Molecule_Maker(
            2,
            name=['CA', 'P'],
            resname=['ALA', 'DA'],
            resid=[1, 2],
            segname='MIXD',
            moltype=['protein', 'dna'])

        report = molecule.moltype_by_segname_report()

        self.assertEqual(report['overall_status'], 'clean')
        segment = report['segments']['MIXD']
        self.assertEqual(segment['status'], 'resolved_dna')
        self.assertEqual(segment['assigned_moltypes'], ['dna'])

    def test_conflicting_dna_and_rna_evidence_is_reported(self):
        molecule = system.Molecule_Maker(
            2,
            name=["O2'", 'P'],
            resname=['ADE', 'DA'],
            resid=[1, 2],
            segname='HYBR',
            moltype=['nucleic', 'nucleic'])

        report = molecule.moltype_by_segname_report()

        self.assertEqual(report['overall_status'], 'conflict')
        segment = report['segments']['HYBR']
        self.assertEqual(segment['status'], 'conflict')
        self.assertEqual(segment['decision'], 'assigned recognized nucleic atoms as nucleic')
        self.assertEqual(segment['rna_atom_evidence'], ["O2'"])
        self.assertEqual(segment['dna_resname_evidence'], ['DA'])

    def test_blank_segname_report_falls_back_to_chain_groups(self):
        molecule = system.Molecule_Maker(
            3,
            name=["O2'", 'P', 'P'],
            resname=['ADE', 'GUA', 'THY'],
            resid=[1, 2, 3],
            chain=['R', 'R', 'D'],
            segname=['', '', ''],
            moltype=['rna', 'rna', 'dna'])

        report = molecule.moltype_by_segname_report()

        self.assertEqual(report['overall_status'], 'ambiguous')
        self.assertEqual(sorted(report['segments'].keys()), ['chain:D', 'chain:R'])
        self.assertEqual(report['segments']['chain:R']['grouping_source'], 'chain')
        self.assertEqual(report['segments']['chain:R']['chain'], 'R')
        self.assertEqual(report['segments']['chain:D']['ambiguous_resnames'], ['THY'])


if __name__ == '__main__':
    unittest.main()
