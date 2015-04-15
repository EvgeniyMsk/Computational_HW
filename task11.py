from sympy import *
from sympy.matrices import Matrix, diag
import itertools

def find_max_not_diag(a):
    n = a.shape[0]

    pairs = filter(lambda (i, j): i != j, itertools.product(range(0, n), repeat=2))
    elems = [a[i, j] for (i, j) in pairs]
    return max(elems)

def calc_d_cos_sin(a, i, j):
    d     = sqrt((a[i, i] - a[j, j])**2 + 4 * a[i,j]**2)
    cosph = sqrt(0.5 * (1 + abs(a[i, i] - a[j, j]) / d))
    sinph = sqrt(sign(a[i, j] * (a[i, i] - a[j, j])) * sqrt(0.5 * (1 - abs(a[i, i] - a[j, j]) / d)))
    return (d, cosph, sinph)

def build_T(a, i, j):
    n = a.shape[0]
    T = diag([1] * n)

    d, cosph, sinph = calc_d_cos_sin(a, i, j)

    T[i, i] = cosph
    T[i, j] = -sinph
    T[j, i] = sinph
    T[j, j] = cosph

    return T

def matrix_jacob(a):
    pass

def build_C(a, i, j):
    n = a.shape[0]
    c = a

    d, cosph, sinph = calc_d_cos_sin(a, i, j)

    for k in range(0, n):
        c[k, i] = c[i, k] =  cosph * a[k, i] + sinph * a[k, j]
        c[k, j] = c[j, k] = -sinph * a[k, i] + cosph * a[k, j]

    c[i, i] = (a[i, i] + a[j, j] + d * sign(a[i, i] - a[j, j])) / 2
    c[j, j] = (a[i, i] + a[j, j] - d * sign(a[i, i] - a[j, j])) / 2

    return c

if __name__ == "__main__":

    a = Matrix([[1.32, 6.45, 0.38],
                [6.45, 1.51, 0.56],
                [0.38, 0.56, 0.72]])




