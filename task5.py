from sympy import *
from sympy.mpmath import e, quad

x = Symbol('x')

f = 1 / (1.28 - e**(-x-0.01))
f_im = lambda t: f.subs(x, t)

seg = [0, 0.4]

js = quad(f_im, seg)

def rect_integral(func, seg, n):
    h = (seg[1] - seg[0]) / n

    intgr = 0
    for k in range(1, n+1):
        intgr += h * func.subs(x, seg[0] + (2*k-1) * h/2)

    fdiff = func.diff(x, 2)
    fdiffm = max(fdiff.subs(x, seg[0]), fdiff.subs(x, seg[1]))
    error = fdiffm * (seg[1] - seg[0])**3 / (24 * n**2)

    return intgr, error

def trap_integral(func, seg, n):
    h = (seg[1] - seg[0]) / n

    intgr = h/2 * (func.subs(x, seg[0]) + func.subs(x, seg[1]))
    for k in range(1, n):
        intgr += h * func.subs(x, seg[0] + k*h)
    #intgr *= h

    fdiff = func.diff(x, 2)
    fdiffm = max(fdiff.subs(x, seg[0]), fdiff.subs(x, seg[1]))
    error = - fdiffm * (seg[1] - seg[0])**3 / (12 * n**2)

    return intgr, error

def simpson_integral(func, seg, n):
    h = (seg[1] - seg[0]) / n

    intgr = h/3 * (func.subs(x, seg[0]) + func.subs(x, seg[1]))
    for k in range(1, n/2 + 1):
        intgr += h/3 * 4 * func.subs(x, seg[0] + (2*k-1)*h)
    for k in range(1, n/2):
        intgr += h/3 * 2 * func.subs(x, seg[0] + (2*k)*h)

    fdiff = func.diff(x, 4)
    fdiffm = max(fdiff.subs(x, seg[0]), fdiff.subs(x, seg[1]))
    error = - fdiffm * (seg[1] - seg[0])**5 / (2880 * n**4)

    return intgr, error

def runge(j8, j16, k):
    return (2**k * j16 - j8) / (2**k - 1)

def print_integrals(intgr_formula, k):

    j8, A = intgr_formula(f, seg, 8)
    j16   = intgr_formula(f, seg, 16)[0]
    jru   = runge(j8, j16, k)
    err   = abs(js - j8)

    print \
        indent + "{:+.8f}".format(float(j8))   + indent + "|" \
        + indent + "{:+.8f}".format(float(j16))  + indent + "|" \
        + indent + "{:+.8f}".format(float(jru))   + indent + "|" \
        + indent + "{:+.8f}".format(float(A)) + indent + "|" \
        + indent + "{:+.8f}".format(float(err)) + indent + "|"

if __name__ == "__main__":

    print "Integral = {}".format(js)
    print ""

    indent = " " * 4
    print \
      indent + "     J8    " + indent + "|" \
    + indent + "    J16    " + indent + "|" \
    + indent + "  J Runge  " + indent + "|" \
    + indent + "     A     " + indent + "|" \
    + indent + "  |J - J8| " + indent + "|"
    print "-" * 100

    print_integrals(rect_integral, 2)
    print_integrals(trap_integral, 2)
    print_integrals(simpson_integral, 4)
