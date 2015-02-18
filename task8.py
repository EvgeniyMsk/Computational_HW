from sympy import *
from sympy.matrices import *

def gauss_step(a, k):
    irow = 0
    for i in range(k, a.shape[0]):
        if a[i, k] != 0:
            irow = i
            break

    a.row_swap(irow, k)

    a[k, :] *= 1 / a[k, k]

    for i in range(k+1, a.shape[0]):
        a[i, :] -= a[i, k] * a[k, :]

def gauss_direct(a):
    for k in range(0, a.shape[0]):
        gauss_step(a, k)

def gauss_reverse(a):
    nrow = a.shape[0]
    ncol = a.shape[1]
    b = a[:, ncol - 1]
    x = zeros(nrow, 1)
    for k in range(nrow - 1, -1, -1):
        x[k] = b[k]
        for j in range(k + 1, nrow):
            x[k] -= x[j] * a[k, j]
    return x

def check_solution(a, x):
    nrow = a.shape[0]
    ncol = a.shape[1]
    b = a[:, ncol - 1]
    a = a[:, 0:ncol - 1]
    ax = a * x
    for (axi, bi) in zip(ax, b):
        if axi != bi:
            return False
    return True

def calc_sum(a):
    nrow = a.shape[0]
    ncol = a.shape[1]
    b = zeros(nrow, 1)
    for i in range(0, nrow):
        for j in range(0, ncol):
            b[i, 0] += a[i, j]
    return b


if __name__ == "__main__":

    a = Matrix([[2, 2, 4, 0], [5, 2, 9, -2], [3, 12, 0, -1]])

    a_orig = a

    nrow = a.shape[0]
    ncol = a.shape[1]

    print "\n Given matrix: \n"
    pprint(a)

    # calculate checksum
    b = calc_sum(a)
    print " \n Checksum: \n"
    pprint(b)
    # insert it to matrix
    a = a.col_insert(ncol, b)
    # do gauss direct pass
    gauss_direct(a)
    # remove checksum column
    b = a[:, ncol]
    a = a[:, 0:ncol]
    print "\n Gauss direct: \n"
    pprint(a)
    print "\n Checksum column: \n"
    pprint(b)
    # calculate new checksum
    b = calc_sum(a)
    print "\n New checksum column: \n"
    pprint(b)
    # calculate x vector
    x = gauss_reverse(a)
    print "\n Solution of linear system: \n"
    pprint(x)
    print "\n Check solution result: " + str(check_solution(a_orig, x))


