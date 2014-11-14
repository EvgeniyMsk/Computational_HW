from matplotlib import pyplot
from sys import stdin
from sympy import *


x = Symbol('x')
f = sin(Rational(1,3) * x)
f_im = lambda xs: [N(f.subs(x, xi)) for xi in xs]
seg = [0, pi]

def get_lagrange(xs, fs):
    w = 1
    for xi in xs:
        w *= (x - xi)
    dw = w.diff(x)
    lagr = 0
    for (xi, fi) in zip(xs, fs):
        lagr += w / ((x - xi) * dw.subs(x, xi)) * fi
    return (w, lagr)

def print_table(xs, fs, ps):
    w, lagr = get_lagrange(xs, fs)

    fps = f_im(ps)
    ls  = [N(lagr.subs(x, p)) for p in ps]
    difs = [abs(a - b) for (a, b) in zip(fps, ls)]

    df = f.diff(x, len(xs))
    factr = factorial(len(xs))
    df_max = max(abs(df.subs(x, seg[0])), abs(df.subs(x, seg[1])))

    errs = [(abs(w.subs(x, p) * df_max) * Rational(1, factr)) for p in ps]

    indent = " " * 4
    print \
      indent + "     x     " + indent + "|" \
    + indent + "     f     " + indent + "|" \
    + indent + "     L     " + indent + "|" \
    + indent + "   L - f   " + indent + "|" \
    + indent + "     A     " + indent + "|"
    print "-" * 100

    for i in range(0, len(ps)):
        print \
          indent + "{:+.8f}".format(float(ps[i]))   + indent + "|" \
        + indent + "{:+.8f}".format(float(fps[i]))  + indent + "|" \
        + indent + "{:+.8f}".format(float(ls[i]))   + indent + "|" \
        + indent + "{:+.8f}".format(float(difs[i])) + indent + "|" \
        + indent + "{:+.8f}".format(float(errs[i])) + indent + "|"

    pyplot.plot(ps, fps)
    pyplot.plot(ps, ls)
    pyplot.legend(['f', 'L'], loc='upper left')
    pyplot.show()


if __name__ == "__main__":

    xs = [None] * 8
    fs = [None] * 8

    xs[0] = [pi/7, pi/5, pi/3, pi/2, 3*pi/4]
    fs[0] = f_im(xs[0])

    xs[1] = [pi/10, pi/6, pi/3, pi/4, 3*pi/8]
    fs[1] = f_im(xs[1])

    xs[2] = [6*pi/10, 2*pi/3, 3*pi/4, 7*pi/8, 9*pi/10]
    fs[2] = f_im(xs[2])

    xs[3] = [pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4]
    fs[3] = f_im(xs[3])

    xs[4] = xs[0] + [pi/11, pi/9, 6*pi/11, 2*pi/3, 12*pi/13]
    xs[4].sort()
    fs[4] = f_im(xs[4])

    xs[5] = xs[1] + [pi/11, pi/7, 3*pi/10, 2*pi/5, 5*pi/12]
    xs[5].sort()
    fs[5] = f_im(xs[5])

    xs[6] = xs[2] + [4*pi/7, 5*pi/9, 5*pi/6, 7*pi/10, 10*pi/11]
    xs[6].sort()
    fs[6] = f_im(xs[6])

    xs[7] = xs[3] + [3*pi/8, 5*pi/12, 4*pi/7, 7*pi/12, 5*pi/8]
    xs[7].sort()
    fs[7] = f_im(xs[7])

    ps = [(i * pi)/20 for i in range(0, 20)]

    for i in range(0, 9):
        print_table(xs[i], fs[i], ps)

        print "press any key to continue"
        line = stdin.readline()




