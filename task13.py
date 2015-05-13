import sys
import os
import math
from sympy import *
from sympy.mpmath import e

x = Symbol('x')
t = Symbol('t')

KMOD = 100
KMAX = 5
STEP = 10

def f_range(a, b, d):
    while a < b + d:
        yield a
        a += d

def print_header(n, h, out):
    if not out:
        return

    out.write('  {:^8}  |'.format('k'))
    out.write('{indent}{num:^15}{indent}|'.format(indent=' '*4, num='eps'))
    for i in range(0, n+1):
        out.write('{indent}{num:+15.2f}{indent}|'.format(indent=' '*4, num=float(i*h)))
    out.write(os.linesep)

def print_row(k, eps, row, out):

    if k == -1:
        n = 10
        print_header(n, float(1) / n, out)
        return

    fmt_str = '{indent}{num:+15.12f}{indent}|'

    out.write('  {:8.6f}  |'.format(float(k)))
    out.write(fmt_str.format(indent=' '*4, num=float(eps)))
    for el in row:
        out.write(fmt_str.format(indent=' '*4, num=float(el)))
    out.write(os.linesep)

def print_solution(u, a, b, N, K, out):

    h = float(b - a) / N
    tau = h**2 / 2

    n = 10
    print_header(n, float(1) / n, out)
    for k in range(0, K):

        if k % KMOD == 0:
            uk = u.subs(t, k * tau)
            row = [uk.subs(x, xi) for xi in f_range(a, b, h)][:N+1]
            print_row(k * tau, 0, row[0::STEP], out)

def update_err(x_range, row, func):
    a, b, h = x_range
    xs = [func.subs(x, xi) for xi in f_range(a, b, h)]
    errs = [abs(z - y) for (z, y) in zip(xs, row)]
    return max(errs)

def explicit_method(u, f, a, b, N, K, out, flg):

    h = float(b - a) / N
    tau = h**2 / 2

    psi_1 = u.subs(x, a)
    psi_2 = u.subs(x, b)
    u0    = u.subs(t, 0)

    coefs = zeros(N+1, N+1)
    for i in range(1, N):
        coefs[i, i-1] = tau / h**2
        coefs[i, i]   = 1 - 2*tau / h**2
        coefs[i, i+1] = tau / h**2

    u_pr = Matrix(N+1, 1, [u0.subs(x, xi) for xi in f_range(a, b, h)][:N+1])

    print_row(-1, 0, u_pr, out)
    if flg:
        print_row(0, 0, u_pr[0::STEP], out)
    else:
        print_row(0, 0, u_pr, out)

    for k in range(1, K):
        fcol = Matrix(N+1, 1, [f.subs([(t, k * tau), (x, xi)]) for xi in f_range(a, b, h)][:N+1])
        fcol *= tau
        fcol[0] = psi_1.subs(t, k * tau)
        fcol[N] = psi_2.subs(t, k * tau)

        u_new = coefs * u_pr + fcol
        err = update_err((a, b, h), u_new, u.subs(t, k * tau))

        if flg:
            if k % KMOD == 0:
                print_row(k * tau, err, u_new[0::STEP], out)
        else:
            if k <= KMAX:
                print_row(k * tau, err, u_new, out)

        u_pr = u_new

        if k % 10 == 0:
            print 'Finish {} rows'.format(k)

def implicit_method(u, f, a, b, N, K, out, flg):

    h = float(b - a) / N
    tau = h**2 / 2

    psi_1 = u.subs(x, a)
    psi_2 = u.subs(x, b)
    u0    = u.subs(t, 0)

    coefs = zeros(N+1, N+1)
    coefs[0, 0] = coefs[N, N] = 1
    for i in range(1, N):
        coefs[i, i-1] = tau / h**2
        coefs[i, i]   = -1 - 2*tau / h**2
        coefs[i, i+1] = tau / h**2

    usym = [Symbol('u' + str(i)) for i in range(0, N+1)]

    u_pr = Matrix(N+1, 1, [u0.subs(x, xi) for xi in f_range(a, b, h)][:N+1])

    print_row(-1, 0, u_pr, out)
    if flg:
        print_row(0, 0, u_pr[0::STEP], out)
    else:
        print_row(0, 0, u_pr, out)

    for k in range(1, K):
        fcol = Matrix(N+1, 1, [(-tau * f - u_pr[i]).subs([(t, k * tau), (x, xi)]) for (i, xi) in zip(range(0, N+1), f_range(a, b, h))])
        fcol[0] = psi_1.subs(t, k * tau)
        fcol[N] = psi_2.subs(t, k * tau)
        coefs = coefs.col_insert(N+1, fcol)

        sol = solve_linear_system(coefs, *usym)
        u_new = [sol[ui] for ui in usym]
        err = update_err((a, b, h), u_new, u.subs(t, k * tau))

        if flg:
            if k % KMOD == 0:
                print_row(k * tau, err, u_new[0::STEP], out)
        else:
            if k <= KMAX:
                print_row(k * tau, err, u_new, out)

        u_pr = u_new

        if k % 10 == 0:
            print 'Finish {} rows'.format(k)

        coefs.col_del(N+1)



if __name__ == "__main__":

    u = e**(-t/4) * sin(x/2) + e**(-t) * (1 - x**2)

    f = e**(-t) * (1 + x**2)

    a = 0
    b = 1
    K = 500

    N = 10

    out = open('explicit_method.txt', 'w')

    print 'Explicit method is performed [N={}]...'.format(N)

    out.write('Explicit method N={N}{sep}'.format(N=N, sep=os.linesep))
    explicit_method(u, f, a, b, N, K, out, False)
    out.write(os.linesep)

    N = 100

    print 'Explicit method is performed [N={}]...'.format(N)

    out.write('Explicit method N={N}{sep}'.format(N=N, sep=os.linesep))
    explicit_method(u, f, a, b, N, K, out, True)
    out.write(os.linesep)

    print "Calculate solution ..."

    out.write('Solution{}'.format(os.linesep))
    print_solution(u, a, b, N, K, out)


    N = 10

    out = open('implicit_method.txt', 'w')

    print 'Implicit method is performed [N={}]...'.format(N)

    out.write('Implicit method N={N}{sep}'.format(N=N, sep=os.linesep))
    implicit_method(u, f, a, b, N, K, out, False)
    out.write(os.linesep)

    N = 100св 

    print 'Implicit method is performed [N={}]...'.format(N)

    out.write('Implicit method N={N}{sep}'.format(N=N, sep=os.linesep))
    implicit_method(u, f, a, b, N, K, out, True)
    out.write(os.linesep)

    print "Calculate solution ..."

    out.write('Solution{}'.format(os.linesep))
    print_solution(u, a, b, N, K, out)



