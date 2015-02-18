import os
from sympy import *
from sys import stdout
from matplotlib import pyplot


x = Symbol('x')
y = Symbol('y')
a = Symbol('a')
b = Symbol('b')

# f = (a*y**2 + x**2 + 2)/(b*y + x + 4)
f = -y * (1+x)

def f_range(a, b, d):
    while a < b:
        yield a
        a += d

# def differences(f, h, x0, y0):
#     hs = [h * f.subs([(x, x0), (y, y0)])]
#     xi = x0
#     yi = y0
#     while True:
#         yield hs
#         xi = xi + h
#         yi = yield
#         hi = h * f.subs([(x, xi), (y, yi)])
#         hsi = [hi]
#         #print hs
#         for hk in hs:
#             hi = hi - hk
#             hsi.append(hi)
#         hs = hsi

def differences(f, h, xi, yi, prev):
    et = h * f.subs([(x, xi), (y, yi)])
    ets = [et]
    for pr in prev:
        et = et - pr
        ets.append(et)
    return ets

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

        print "{:+.8f}".format(float(k1))
        print "{:+.8f}".format(float(k2))
        print "{:+.8f}".format(float(k3))
        print "{:+.8f}".format(float(k4))
        print ""

        yk_1 = yk + Rational(1, 6) * (k1 + 2*k2 + 2*k3 + k4)
        ps.append((xk + h, yk_1))
        yk = yk_1
    return ps[:n+1]

def adams_method(f, x0, y0, xn, n):
    h = float(xn - x0) / n
    ps4 = runge_kutt_method(f, x0, y0, x0 + 4*h, 4)

    diffs = []
    diff = []
    for p in ps4:
        diff = differences(f, h, p[0], p[1], diff)
        diffs.append(diff)


    indent = " " * 4
    print \
      indent + "   xk  " + indent + "|" \
    + indent + "   yk  " + indent + "|" \
    + indent + "   h   " + indent + "|" \
    + indent + "  d1h  " + indent + "|" \
    + indent + "  d2h  " + indent + "|" \
    + indent + "  d3h  " + indent + "|" \
    + indent + "  d4h  " + indent + "|"
    print "-" * 112

    for (p,diff)  in zip(ps4, diffs):
        line = list(p) + diff
        for diffi in line:
            stdout.write("  {:+.8f}  |".format(float(diffi)))
        stdout.write(os.linesep)

    ps = ps4
    xi = ps4[4][0]
    yi = ps4[4][1]
    idiff = diffs[4]
    for i in range(5, 11):
        xi += h
        if xi > xn:
            break

        yi = yi + idiff[0] + Rational(1, 2) * idiff[1] + Rational(5, 12) * idiff[2] + \
            Rational(3, 8) * idiff[3] + Rational(251, 720) * idiff[4]

        #print "y{i} = {yi:+.6f} + {diff0:+.6f} + 1/2 * {diff1:+.6f} + 5/12 * {diff2:+.6f} + 3/8 * {diff3:+.6f} + 251/720 * {diff4:+.6f}".format(
        #    i=i, yi=yi, diff0=idiff[0], diff1=idiff[1], diff2=idiff[2], diff3=idiff[3], diff4=idiff[4]
        #)
        #i += 1

        ps.append((xi, yi))

        idiff = differences(f, h, xi, yi, idiff)
        line = [xi, yi] + idiff
        for (k, diff) in zip(range(0, 7), line):
            stdout.write("  {:+.8f}  |".format(float(diff)))
        stdout.write(os.linesep)

    print ""

    return ps

def calc_err(f, ps, h):
    # M coeffs
    M1 = Rational(3, 5)
    M2 = Rational(7, 25)
    M3 = Rational(5, 4)
    M4 = M2 + M1 * M3

    AS = []
    for (xk, yk) in ps:
        A = (M4 / 2*M3) * h * f.subs([(x, xk), (y, yk)]) * exp(M3 * (x0 - xk))
        AS.append(A)
    return AS


def print_table(points, y_pr, err):
    indent = " " * 4
    print \
      indent +     " k "    + indent + "|" \
    + indent + " xk  "       + indent + "|" \
    + indent + "    yk     " + indent + "|" \
    + indent + "  yk prec  " + indent + "|" \
    + indent + "yk - yk_pr " + indent + "|" \
    + indent + "     A     " + indent + "|"
    print "-" * 106

    n = len(points)
    for (k, (xk, yk), errk) in zip(range(0, n), points, err):
        print \
        indent   + "{:<3}".format(k)             + indent + "|"     \
        + indent + "{:+.2f}".format(float(xk))   + indent + "|"     \
        + indent + "{:+.8f}".format(float(yk))   + indent + "|"     \
        + indent + "{:+.8f}".format(float(y_pr(xk))) + indent + "|" \
        + indent + "{:+.8f}".format(float(abs(yk - y_pr(xk))))   + indent + "|"     \
        + indent + "{:+.8f}".format(float(errk)) + indent + "|"

    print "-" * 106


if __name__ == "__main__":

    # f = f.subs([(a, 1), (b, 1)])

    x0 = 0
    y0 = 1
    xn = 1
    n = 10

    g = Function('g')
    y_pr = mpmath.odefun(lambda xi, yi: f.subs([(x, xi), (y, yi)]), 0, 0)

    print "Euler method h = 0.1"
    print ""
    euler_res1 = euler_method(f, x0, y0, xn, n)
    err = calc_err(f, euler_res1, 0.1)
    print_table(euler_res1, y_pr, err)
    print ""

    print "Euler method h = 0.05"
    print ""
    euler_res2 = euler_method(f, x0, y0, xn, 2*n)
    err = calc_err(f, euler_res2, 0.05)
    print_table(euler_res2, y_pr, err)
    print ""

    print "Euler method h = 0.2"
    print ""
    euler_res3 = euler_method(f, x0, y0, xn, n/2)
    err = calc_err(f, euler_res3, 0.2)
    print_table(euler_res3, y_pr, err)
    print ""

    print "Runge-Kutta method"
    print ""
    runge_kutt_res = runge_kutt_method(f, x0, y0, xn, n)
    err = calc_err(f, runge_kutt_res, 0.1)
    print_table(runge_kutt_res, y_pr, err)
    print ""

    print "Adams method"
    print ""
    adams_res = adams_method(f, x0, y0, xn, n)
    err = calc_err(f, adams_res, 0.1)
    print_table(adams_res, y_pr, err)


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