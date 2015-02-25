from sympy import *
from sympy.matrices import *
from cmath import *

from task8 import gauss_reverse


def triangulate_matrix(a):
    u = zeros(a.shape[0], a.shape[1])
    nrow = a.shape[0]
    ncol = a.shape[1]

    for i in range(0, nrow):
        usum = 0
        for k in range(0, i):
            usum += u[k, i]**2
        u[i, i] = sqrt(a[i, i] - usum)
        for j in range(i+1, ncol):
            usum = 0
            for k in range(0, i):
                usum += u[k,i] * u[k,j]
            u[i, j] = (a[i, j] - usum) / u[i,i]

    return u

if __name__ == "__main__":

    a = Matrix([[1.25, 0.34, 0.01, -1.27],
                [0.34, -1.68, -2.04, 0.25],
                [0.01, -2.04, 0.94, 2.45],
                [-1.27, 0.25, 2.45, 0.85]])

    b = Matrix([[4.241942,
                 -0.030203,
                 -7.411342,
                 -6.289385]]).T

    nrow = a.shape[0]
    ncol = a.shape[1]

    print "\n Given matrix A: \n"
    pprint(a)

    print "\n Given column b: \n"
    pprint(b)

    u = triangulate_matrix(a)
    print "\n Matrix U: \n"
    pprint(u)

    ut = transpose(u)
    print "\n U^T * U: \n"
    pprint(ut * u)

    ut = ut.col_insert(ncol, b)
    y = gauss_reverse(ut, False)
    print "\n y: \n"
    pprint(y)

    u = u.col_insert(ncol, y)
    x = gauss_reverse(u)
    print "\n x: \n"
    pprint(x)

    print "\n Ax: \n"
    pprint(a * x)

    print "\n Given column b: \n"
    pprint(b)