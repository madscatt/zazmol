from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

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
#	LINEAR_ALGEBRA
#
#	12/13/2009	--	initial coding			:	jc
#	12/24/2015	--	refactored for release  :	jc
#
#	 1         2         3         4         5         6         7
# LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	LINEAR_ALGEBRA contains methods to perform basic mathematical operations
'''

import sys
import numpy
import math
import sasmol.matrix_math as matrix_math

def cross_product(a, b):
    ''' 
    Returns cross product (vector product) on two vectors

    Parameters
        ----------
    a
        float list : vector a

    b
        float list : vector b

    Returns
    -------
    numpy.array of cross product

    Examples
    -------

    >>> import sasmol.linear_algebra as linear_algebra
    >>> a = [1.0, 2.0, 3.0]
    >>> b = [-1.0, 6.0, 8.0]
    >>> linear_algebra.cross_product(a, b)
    array([ -2., -11.,   8.]) 

    '''

    cross = [0] * 3
    cross[0] = a[1] * b[2] - a[2] * b[1]
    cross[1] = a[2] * b[0] - a[0] * b[2]
    cross[2] = a[0] * b[1] - a[1] * b[0]

    return numpy.array(cross)


def matrix_multiply(a, b):
    ''' 
    Returns the result of multiplying matrix a by matrix b

    Parameters
        ----------
    a
        float list : matrix a

    b
        float list : matrix b

    Returns
    -------
    tuple 
        error
            list with error code (if error occurs)

        numpy.array of matrix product

    Examples
    -------

    >>> import sasmol.linear_algebra as linear_algebra
    >>> import numpy
    >>> a=numpy.array([[5.0, 3.0, 1.0],[2.0, 3.0, 5.0]])
    >>> b=numpy.array([[2.0, -4.0, 8.0]]).T
    >>> linear_algebra.matrix_multiply(a, b)
    ([], array([[  6.],
       [ 32.]]))


    '''
    error = []

    shape_a = a.shape
    shape_b = b.shape
    dim_a1 = a.shape[0]
    dim_a2 = a.shape[1]
    dim_b1 = b.shape[0]

    try:
        dim_b2 = b.shape[1]
    except:
        dim_b2 = 1

    c = numpy.zeros((dim_a1, dim_b2), numpy.float)
    if(dim_a2 != dim_b1):
        message = 'incompatible matrices'
        error.append(message)
        return error, c

    c = matrix_math.matrix_multiply(a, b, dim_a1, dim_a2, dim_b2)

    return error, c


def find_u(x, y):
    '''
    Method to find the U matrix used to align two sets of 3 x 3 coordinate
    arrays

    Parameters
        ----------
    x
        numpy array : 3 x 3

    y
        numpy array : 3 x 3

    Returns
    -------
    numpy.array
        containing U matrix for alignment of two vectors

    Examples
    -------

    >>> import sasmol.linear_algebra as linear_algebra
    >>> import numpy
    >>> a = numpy.array([[[1.0, 2.0, 3.0],[1.0, 2.0, 3.0],[1.0, 2.0, 3.0]],[[1.0, 2.0, 3.0],[1.0, 2.0, 3.0],[1.0, 2.0, 3.0]],[[1.0, 2.0, 3.0],[1.0, 2.0, 3.0],[1.0, 2.0, 3.0]]])
    >>> b = numpy.array([[[-1.0, 6.0, 8.0],[-1.0, 6.0, 8.0],[-1.0, 6.0, 8.0]],[[-1.0, 6.0, 8.0],[-1.0, 6.0, 8.0],[-1.0, 6.0, 8.0]],[[-1.0, 6.0, 8.0],[-1.0, 6.0, 8.0],[-1.0, 6.0, 8.0]]])
    >>> b = [-1.0, 6.0, 8.0]
    >>> linear_algebra.find_u(a, b)
    array([[-0.33333333, -0.33333333, -0.33333333],
        [-0.33333333, -0.33333333, -0.33333333],
        [-0.33333333, -0.33333333, -0.33333333]])

    '''

    b = numpy.zeros(9, numpy.float)
    k = 0
    for i in range(3):
        for j in range(3):
            rad = 0.0
            for n in range(len(x)):
                rad = rad + y[n, i] * x[n, j]
            numpy.put(b, k, rad)
            k = k + 1
    r = numpy.reshape(b, (-1, 3))
    r = numpy.mat(r)
    rt = r.T		# transpose of r
    rtr = rt * r  # matrix multiply rt * r
    uk, ak = numpy.linalg.eig(rtr)
    ak = ak.T

    idx = uk.argsort()
    idx = idx[::-1]
    uk = uk[idx]
    ak = ak[idx]

    ak[2] = numpy.cross(ak[0], ak[1])
    rak0 = numpy.inner(r, ak[0])
    rak1 = numpy.inner(r, ak[1])
    rak0.shape = (1, 3)
    rak1.shape = (1, 3)

    if(uk[0] == 0.0):
        urak0 = (10**15) * rak0
    else:
        urak0 = (1.0 / math.sqrt(abs(uk[0]))) * rak0

    if(uk[1] == 0.0):
        urak1 = (10**15) * rak1
    else:
        urak1 = (1.0 / math.sqrt(abs(uk[1]))) * rak1

    urak2 = numpy.cross(urak0, urak1)
    bk = numpy.zeros((3, 3), numpy.float)
    bk[0] = urak0
    bk[1] = urak1
    bk[2] = urak2
    lu = numpy.zeros(9, numpy.float)
    lk = 0
    for j in range(3):
        for i in range(3):
            rad = 0.0
            for k in range(3):
                rad = rad + bk[k, i] * ak[k, j]
            numpy.put(lu, lk, rad)
            lk = lk + 1
    u = numpy.reshape(lu, (-1, 3))

    return u


def signed_angle(a, b, c):
    '''
    This method calcultes the sign of the angle which is used in the calculation of a dihedral angle.
    As such, this method is not validated for other uses.  It will fail if the basis atoms for the
    dihedral (atom 2 and atom 3) overlap.

    Parameters
        ----------
    a
        float list : vector a

    b
        float list : vector b

    c
        float list : vector c

    Returns
    -------
    float
        signed angle in degrees

    Examples
    -------

    >>> import sasmol.linear_algebra as linear_algebra
    >>> a = [1.0, 2.0, 3.0]
    >>> b = [-1.0, 6.0, 8.0]
    >>> c = [-4.0, -1.0, 4]
    >>> linear_algebra.signed_angle(a, b, c)
    21.444512921997863   

    '''

    ada = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    bdb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2])

    if(ada * bdb <= 0.0):
        return 180.0
    else:
        adb = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
        try:
            argument = adb / math.sqrt(ada * bdb)
            angle = (180.0 / math.pi) * math.acos(argument)
        except:
            return 180.0

    cp = cross_product(a, b)
    dp = cp[0] * c[0] + cp[1] * c[1] + cp[2] * c[2]
    sign = cmp(dp, 0.0)

    return sign * angle


def dihedral_angle(a1, a2, a3, a4):
    '''
    Calculates the dihedral angle between four vectors  


    Parameters
        ----------
    a1
        float list : vector a1

    a2
        float list : vector a2

    a3
        float list : vector a3

    a4
        float list : vector a4

    Returns
    -------
    float
        dihedral angle in degrees

    Examples
    -------

    >>> import sasmol.linear_algebra as linear_algebra
    >>> a1 = numpy.array([1.0, 2.0, 3.0])
    >>> a2 = numpy.array([-1.0, 6.0, 8.0])
    >>> a3 = numpy.array([-4.0, -1.0, 4.0])
    >>> a4 = numpy.array([-3.0, -41, 3.0])
    >>> linear_algebra.dihedral_angle(a1, a2, a3, a4)
    85.950635659264 

    '''

    r1 = numpy.zeros(3, numpy.float)
    r2 = numpy.zeros(3, numpy.float)
    r3 = numpy.zeros(3, numpy.float)
    r4 = numpy.zeros(3, numpy.float)

    r1 = a1 - a2
    r2 = a3 - a2
    r3 = -1.0 * r2
    r4 = a4 - a3

    n1 = cross_product(r1, r2)
    n2 = cross_product(r3, r4)

    dihedral_angle = signed_angle(n1, n2, r2)

    return dihedral_angle


def calculate_angle(a, b, c):
    '''

    Calculates the dihedral angle between three vectors  

    Parameters
        ----------
    a
        float list : vector a

    b
        float list : vector b

    c
        float list : vector c

    Returns
    -------
    float
        angle in radians

    Examples
    -------

    >>> import sasmol.linear_algebra as linear_algebra
    >>> a = numpy.array([1.0, 2.0, 3.0])
    >>> b = numpy.array([-1.0, 6.0, 8.0])
    >>> c = numpy.array([-4.0, -1.0, 4.0])
    >>> linear_algebra.calculate_angle(a, b, c)
    0.7556508878558726 

    '''
    u = a - b
    v = c - b

    norm_u = math.sqrt(sum(u * u))
    norm_v = math.sqrt(sum(v * v))
    c = numpy.dot(u, v) / (norm_u * norm_v)  # -> cosine of the angle
    angle = math.acos(c)

    return angle
