import os
from sympy import *
from sys import stdout
from matplotlib import pyplot


x = Symbol('x')
y = Symbol('y')
a = Symbol('a')
b = Symbol('b')

f = (a*y**2 + x**2 + 2)/(b*y + x + 4)

def f_range(a, b, d):
    while a < b:
        yield a
        a += d

def differences(f, h, x0, y0):
    hs = [h * f.subs([(x, x0), (y, y0)])]
    xi = x0
    yi = y0
    while True:
        yield hs
        xi = xi + h
        yi = yield
        hi = h * f.subs([(x, xi), (y, yi)])
        hsi = [hi]
        #print hs
        for hk in hs:
            hi = hi - hk
            hsi.append(hi)
        hs = hsi


def euler_method(f, x0, y0, xn, n):
    h = float(xn - x0) / n
    ps = [(x0, y0)]
    yk = y0
    for xk in f_range(x0, xn, h):
        yk_1 = yk + h * f.subs([(x, xk), (y, yk)])
        ps.append((xk + h, yk_1))
        yk = yk_1
    return ps[:n+1]

def runge_kutt_method(f, x0, y0, xn, n):
    h = float(xn - x0) / n
    ps = [(x0, y0)]
    yk = y0
    for xk in f_range(x0, xn, h):
        k1 = h*f.subs([(x, xk), (y, yk)])
        k2 = h*f.subs([(x, xk + h/2), (y, yk + k1/2)])
        k3 = h*f.subs([(x, xk + h/2), (y, yk + k2/2)])
        k4 = h*f.subs([(x, xk + h), (y, yk + k3)])

        yk_1 = yk + Rational(1, 6) * (k1 + 2*k2 + 2*k3 + k4)
        ps.append((xk + h, yk_1))
        yk = yk_1
    return ps[:n+1]

def adams_method(f, x0, y0, xn, n):
    h = float(xn - x0) / n
    ps4 = runge_kutt_method(f, x0, y0, x0 + 4*h, 4)
    gen = differences(f, h, x0, y0)
    diffs = gen.next()
    idiff = []


    indent = " " * 4
    print \
      indent + "   xk  " + indent + "|" \
    + indent + "   yk  " + indent + "|" \
    + indent + "   h   " + indent + "|" \
    + indent + "  d1h  " + indent + "|" \
    + indent + "  d2h  " + indent + "|" \
    + indent + "  d3h  " + indent + "|" \
    + indent + "  d4h  " + indent + "|"
    print "-" * 120

    for p in ps4:
        line = list(p) + diffs
        for diff in line:
            stdout.write("  {:+.8f}  |".format(float(diff)))
        stdout.write(os.linesep)
        gen.next()
        diffs = gen.send(p[1])
        idiff = diffs

    ps = ps4
    xi = ps4[4][0]
    yi = ps4[4][1]
    for none in gen:
        xi += h
        if xi > xn:
            break
        yi = yi + idiff[0] + Rational(1, 2) * idiff[1] + Rational(5, 12) * idiff[2] + \
            Rational(3, 8) * idiff[3] + Rational(251, 720) * idiff[4]
        ps.append((xi, yi))

        idiff = gen.send(yi)
        line = [xi, yi] + idiff
        for (k, diff) in zip(range(0, 7), line):
            stdout.write("  {:+.8f}  |".format(float(diff)))
        stdout.write(os.linesep)

    print ""

    return ps



def print_table(points, y_pr):
    indent = " " * 4
    print \
      indent +     " k "    + indent + "|" \
    + indent + " xk  "       + indent + "|" \
    + indent + "    yk     " + indent + "|" \
    + indent + "  yk prec  " + indent + "|"
    print "-" * 100

    n = len(points)
    for (k, (xk, yk)) in zip(range(0, n), points):
        print \
        indent   + "{:<3}".format(k)             + indent + "|"     \
        + indent + "{:+.2f}".format(float(xk))   + indent + "|"     \
        + indent + "{:+.8f}".format(float(yk))   + indent + "|"     \
        + indent + "{:+.8f}".format(float(y_pr(xk))) + indent + "|"

    print "-" * 100


if __name__ == "__main__":

    f = f.subs([(a, 1), (b, 10)])

    x0 = 0
    y0 = 0
    xn = 1
    n = 10

    g = Function('g')
    y_pr = mpmath.odefun(lambda xi, yi: f.subs([(x, xi), (y, yi)]), 0, 0)

    print "Euler method h = 0.1"
    print ""
    euler_res1 = euler_method(f, x0, y0, xn, n)
    print_table(euler_res1, y_pr)
    print ""

    print "Euler method h = 0.05"
    print ""
    euler_res2 = euler_method(f, x0, y0, xn, 2*n)
    print_table(euler_res2, y_pr)
    print ""

    print "Euler method h = 0.2"
    print ""
    euler_res3 = euler_method(f, x0, y0, xn, n/2)
    print_table(euler_res3, y_pr)
    print ""

    print "Runge-Kutta method"
    print ""
    runge_kutt_res  = runge_kutt_method(f, x0, y0, xn, n)
    print_table(runge_kutt_res, y_pr)
    print ""

    print "Adams method"
    print ""
    adams_res       = adams_method(f, x0, y0, xn, n)
    print_table(adams_res, y_pr)


    xs = [x for (x, y) in euler_res1]
    ys = [y for (x, y) in euler_res1]
    pyplot.plot(xs, ys)

    xs = [x for (x, y) in euler_res2]
    ys = [y for (x, y) in euler_res2]
    pyplot.plot(xs, ys)

    xs = [x for (x, y) in euler_res3]
    ys = [y for (x, y) in euler_res3]
    pyplot.plot(xs, ys)

    xs = [x for (x, y) in runge_kutt_res]
    ys = [y for (x, y) in runge_kutt_res]
    pyplot.plot(xs, ys)

    xs = [x for (x, y) in adams_res]
    ys = [y for (x, y) in adams_res]
    pyplot.plot(xs, ys)

    xs = [x for (x, y) in adams_res]
    ys = [y_pr(x) for x in xs]
    pyplot.plot(xs, ys)

    pyplot.legend(['Euler 0.1', 'Euler 0.05', 'Euler 0.2', 'Runge-Kutta', 'Adams', "Y"], loc='upper left')
    pyplot.show()