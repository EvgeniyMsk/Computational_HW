from sympy import *


x = Symbol('x')
y = Symbol('y')
a = Symbol('a')
b = Symbol('b')

f = (a*y**2 + x**2 + 2)/(b*y + x + 4)

def f_range(a, b, d):
    while a < b:
        yield a
        a += d

def euler_method(f, x0, y0, xn, n):
    h = float(xn - x0) / n
    ps = [(x0, y0)]
    yk = y0
    for xk in f_range(x0, xn, h):
        yk_1 = yk + h * f.subs([(x, xk), (y, yk)])
        ps.append((xk + h, yk_1))
        yk = yk_1
    return ps

def runge_cook_method(f, x0, y0, xn, n):
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
    return ps

def print_table(points, y_pr):
    indent = " " * 4
    print \
      indent +     " k "    + indent + "|" \
    + indent + "    xk     " + indent + "|" \
    + indent + "    yk     " + indent + "|" \
    + indent + "  yk prec  " + indent + "|"
    print "-" * 100

    for (k, (xk, yk)) in zip(range(0, n + 1), points):
        print \
        indent   + "{:<3}".format(k)             + indent + "|"     \
        + indent + "{:+.8f}".format(float(xk))   + indent + "|"     \
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

    euler_res       = euler_method(f, x0, y0, xn, n)
    runge_cook_res  = runge_cook_method(f, x0, y0, xn, n)

    print_table(euler_res, y_pr)
    print_table(runge_cook_res, y_pr)