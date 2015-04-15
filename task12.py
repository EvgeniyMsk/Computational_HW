
from sympy import *
from sympy.mpmath import e

x = Symbol('x')

def test_model(xs, yx, us, vs, kappa, nu):

    space = ' '*5
    print space + ' x ' + space + "|" + space + 'y(x)' + space + "|" + space + 'y(u, v)' + space
    for (i, u, v) in zip(range(0, len(xs)), us, vs):
        yi  = yx.subs(x, xs[i])
        yi_ = u * yx.subs(x, xs[i+1]) + v
        print "{:+12.8f} | {:+12.8f} | {:+12.8f}".format(xs[i], float(yi), float(yi_))
    yn  = yx.subs(x, xs[-1])
    yn_ = kappa * yx.subs(x, xs[-2]) + nu
    print "{:+12.8f} | {:+12.8f} | {:+12.8f}".format(xs[-1], float(yn), float(yn_))

if __name__ == "__main__":

    # solve differential equasition:
    #
    # y'' + p(x)y' + q(x)y = f(x)
    #
    # y'(0) = alpha * y(0)
    # y'(1) = -beta * y(1)

    x0 = 0
    xn = 1
    n = 10
    h = float(xn - x0) / n

    xs = [x0 + i * h for i in range(0, n + 1)]
    ysym = [Symbol('y' + str(i)) for i in range(0, n + 1)]

    alpha = 0.3
    beta  = 1.2

    px = e**x
    qx = x**2 - x - 1
    fx = 1 / (1 + x)

    # model solution, for debugging
    #yx   = sin(x) + 1
    #dyx  = diff(yx, x)
    #d2yx = diff(yx, x, 2)
    #alpha = (dyx / yx).subs(x, 0)
    #beta  = (dyx / yx).subs(x, 1)
    #fx   = d2yx + px * dyx + qx * yx

    a = 1 + px * h/2
    b = 2 - qx * h**2
    c = 1 - px * h/2
    g = fx * h**2

    kappa1 = ((b - 4*a) / (c - 2 * a * h * alpha - 3 * a)).subs(x, xs[1])
    nu1    = (g / (c - 2 * a * h * alpha - 3 * a)).subs(x, xs[1])

    kappa2 = ((4*c - b) / (3 * c - 2 * c * h * beta - a)).subs(x, xs[n-1])
    nu2    = (-g / (3 * c - 2 * c * h * beta - a)).subs(x, xs[n-1])

    us = [0] * (n)
    vs = [0] * (n)
    us[0] = kappa1
    vs[0] = nu1
    for i in range(1, n):
        us[i] = (a / (b - c * us[i-1])).subs(x, xs[i])
        vs[i] = ((c * vs[i-1] - g) / (b - c * us[i-1])).subs(x, xs[i])

    #test_model(xs, yx, us, vs, kappa2, nu2)

    ys = [0] * (n+1)
    ys[n] = float((nu2 + kappa2 * vs[n-1]) / (1 - kappa2 * us[n-1]))

    for i in range(n-1, -1, -1):
        ys[i] = solve(ysym[i] - us[i] * ys[i+1] - vs[i], ysym[i])[0]

    space = ' '*4
    print space + ' x' + space + "|" + space + '  y'
    print '-' * 40
    for (xi, yi) in zip(xs, ys):
        print " {:+6.5f} | {:+12.8f}".format(xi, yi)
