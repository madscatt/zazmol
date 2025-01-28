import unittest
import numpy
#from sasmol.linear_algebra import matrix_multiply
from matrix_math import matrix_multiply  # Import the Fortran extension method

class TestMatrixMultiply(unittest.TestCase):


    def test_matrix_multiply(self):
        # Define test matrices
        a = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=numpy.float32)
        b = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=numpy.float32).T

        # Expected result of the multiplication
        expected_result = numpy.array([[1.0, 0.0], [0.0, 1.0]], dtype=numpy.float32)

        # Perform matrix multiplication
        error, result = matrix_multiply(a, b)
        print('Matrix a:')
        print(a)
        print('Matrix b:')
        print(b)
        print('error', error)
        print('result', result)
        # Check if the result matches the expected result
        numpy.testing.assert_array_almost_equal(result, expected_result)

    '''
    def test_matrix_multiply(self):
        # Define test matrices
        a = numpy.array([[1, 2, 3], [4, 5, 6]], dtype=numpy.float32)
        b = numpy.array([[7, 8], [9, 10], [11, 12]], dtype=numpy.float32)

        a = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=numpy.float32)
        b = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=numpy.float32).T

        # Expected result of the multiplication
        expected_result = numpy.array([[58, 64], [139, 154]], dtype=numpy.float32)
        expected_result = numpy.array([[1.0, 0.0], [0.0, 1.0]], dtype=numpy.float32)

        # Perform matrix multiplication
        error, result = matrix_multiply(a, b)
        print('error', error)
        print('result', result)
        # Check if the result matches the expected result
        numpy.testing.assert_array_almost_equal(result, expected_result)
    '''

if __name__ == '__main__':
    unittest.main()