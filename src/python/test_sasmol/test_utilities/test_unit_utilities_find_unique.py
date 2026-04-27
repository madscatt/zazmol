from unittest import main
import unittest
import sasmol.utilities as utilities

class Test(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(sorted(utilities.find_unique([1,2,2,3])), [1,2,3])

if __name__ == '__main__':
    main()
