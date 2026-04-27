'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D.
'''

import unittest

import sasmol.system as system
import sasmol.utilities as utilities


class Test_unit_utilities_parse_fasta_and_integrity(unittest.TestCase):

    def test_parse_fasta_multi_record_ignores_digits_spaces_and_asterisk(self):
        fasta_lines = [
            '>seq_1\n',
            'A 1 C*G\n',
            ';comment\n',
            'T 2 A*\n',
        ]

        result = utilities.parse_fasta(fasta_lines)

        self.assertEqual(result, ['ACG', 'TA'])

    def test_parse_fasta_returns_error_on_empty_line(self):
        fasta_lines = [
            '>seq_1\n',
            'ACG\n',
            '\n',
            'TTT\n',
        ]

        result = utilities.parse_fasta(fasta_lines)

        self.assertEqual(
            result, 'ERROR: empty lines in fasta sequence are not allowed')

    def test_check_integrity_reports_expected_lengths_for_valid_object(self):
        molecule = system.Molecule_Maker(2)

        result = utilities.check_integrity(molecule, warn=False)

        self.assertEqual(result['_atom'], 2)
        self.assertEqual(result['_index'], 2)
        self.assertEqual(result['_name'], 2)
        self.assertEqual(result['_resid'], 2)

    def test_check_integrity_reports_mismatched_length_in_result_map(self):
        molecule = system.Molecule_Maker(2)
        molecule.setName(['N'])

        result = utilities.check_integrity(molecule, warn=False)

        self.assertEqual(result['_name'], 1)
        self.assertEqual(result['_index'], 2)
        self.assertEqual(result['_resid'], 2)


if __name__ == '__main__':
    unittest.main()
