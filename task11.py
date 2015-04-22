from sympy import *
from sympy.matrices import Matrix, diag
import itertools
import operator

eps = 1e-5

def print_deg_iter(k, ev):

    print "Iteration {k}; Eigen value: {ev:+12.8f}".format(k=k, ev=float(ev))

def deg_method(a, print_cb):
    n = a.shape[0]

    k = 1
    x = ones(n, 1)
    xk  = a * x
    ev = 1e10
    while True:
        xk1 = a * xk
        i = next(i for i in range(0, n) if (xk[i] > eps and xk1[i] > eps))
        ev1 = xk1[i] / xk[i]
        print_cb(k, ev1)

        if (abs(ev1 - ev) < eps):
            break

        ev = ev1
        k += 1
        xk = xk1


def find_max_not_diag(a):
    n = a.shape[0]

    pairs = filter(lambda (i, j): i != j, itertools.product(range(0, n), repeat=2))
    elems = [(abs(a[i, j]), i, j) for (i, j) in pairs]
    return max(elems, key=operator.itemgetter(0))

def calc_d_cos_sin(a, i, j):
    d     = sqrt((a[i, i] - a[j, j])**2 + 4 * a[i,j]**2)
    cosph = sqrt(0.5 * (1 + abs(a[i, i] - a[j, j]) / d))
    sinph = sign(a[i, j] * (a[i, i] - a[j, j])) * sqrt(0.5 * (1 - abs(a[i, i] - a[j, j]) / d))
    return (d, cosph, sinph)

def build_T(a, i, j):
    n = a.shape[0]
    T = eye(n)

    d, cosph, sinph = calc_d_cos_sin(a, i, j)

    T[i, i] = cosph
    T[i, j] = -sinph
    T[j, i] = sinph
    T[j, j] = cosph

    return T

def build_ATA(a, i , j):
    T = build_T(a, i, j)
    return T.T * a * T


def build_C(a, i, j):
    n = a.shape[0]
    c = a.copy()

    d, cosph, sinph = calc_d_cos_sin(a, i, j)

    for k in range(0, n):
        c[k, i] = c[i, k] =  cosph * a[k, i] + sinph * a[k, j]
        c[k, j] = c[j, k] = -sinph * a[k, i] + cosph * a[k, j]

    c[i, i] = (a[i, i] + a[j, j] + d * sign(a[i, i] - a[j, j])) / 2
    c[j, j] = (a[i, i] + a[j, j] - d * sign(a[i, i] - a[j, j])) / 2
    c[i, j] = c[j, i] = 0

    return c

def matrix_jacob(a_, build_cb, print_cb):

    a = a_
    k = 1
    while True:
        max_tuple = find_max_not_diag(a)
        if  max_tuple[0] < eps:
            break
        a = build_cb(a, max_tuple[1], max_tuple[2])
        print_cb(a, k)
        k += 1

def print_jacob_iter(a, k):

    print "    Iteration " + str(k) + '\n'
    pprint(a)
    print '-'*80

if __name__ == "__main__":

    a = Matrix([[1.32, 6.45, 0.38],
                [6.45, 1.51, 0.56],
                [0.38, 0.56, 0.72]])

    print "Degree method (max eigen value) \n"
    deg_method(a, print_deg_iter)

    print '-'*80
    print "\nJacob method (by matrix multiplication) \n"
    matrix_jacob(a, build_ATA, print_jacob_iter)

    print "\nJacob method (by elements modification) \n"
    matrix_jacob(a, build_C, print_jacob_iter)



