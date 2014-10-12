from sympy import *

def f_range(a, b, d):
    while a < b:
        yield a
        a += d

x = Symbol('x')
y = Symbol('y')
a = Symbol('a')
k = Symbol('k')

def filter_near(ps, eps = 0.1):
    res = list()
    l = len(ps)
    for i in range(0, l):
        b = True
        for j in range(i + 1, l):
            pi = ps[i]
            pj = ps[j]
            if (pi[0] - pj[0]) ** 2 + (pi[1] - pj[1]) ** 2 < eps:
                b = False
        if b:
            res.append(ps[i])
    return res

def find_init_approx(f, g, dx, dy, aj, kj, d = 0.05, val = 0.1):
    x0 = dx[0]
    x1 = dx[1]
    y0 = dy[0]
    y1 = dy[1]
    points = list()
    for xi in f_range(x0, x1, d):
        for yi in f_range(y0, y1, d):
            fi = abs(f.subs([(x, xi), (y, yi), (a, aj)]))
            gi = abs(g.subs([(x, xi), (y, yi), (k, kj)]))
            if fi < val and gi < val:
                points.append((xi, yi))
    return filter_near(points)

def newtone_step(f, g, xi, yi):
    dfdx = f.diff(x)
    dfdy = f.diff(y)
    dgdx = g.diff(x)
    dgdy = g.diff(y)
    sub = [(x, xi), (y, yi)]
    mf = -f.subs(sub)
    mg = -g.subs(sub)
    fx =  dfdx.subs(sub)
    fy =  dfdy.subs(sub)
    gx =  dgdx.subs(sub)
    gy =  dgdy.subs(sub)
    detx = matrices.Matrix([mf, fy], [mg, gy]).det()
    dety = matrices.Matrix([fx, mf], [gx, mg]).det()
    det  = matrices.Matrix([fx, fy], [gx, gy]).det()

    xj = xi + float(detx) / det
    yj = yi + float(dety) / det
    return (xj, yj)

def newtone_solve(f, g):
    pass

f = cos(x**2 + y**2) - x + y - a
g = (x + y - 2)**2 / k + (x - y)**2 / (k - 0.1) - 1

print find_init_approx(f, g, (0, 2), (0, 2), 0, 0.4)

#print filter_near([(0, 1), (1, 1), (0, 7)], 2)