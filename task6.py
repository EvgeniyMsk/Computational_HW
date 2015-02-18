from sympy import *
from itertools import combinations

x = Symbol('x')

qs = [1, x, x**2, x**3]

def get_w(n):
    xs = []
    w  = 1
    for i in range(0, n):
        xi = Symbol('x' + str(i))
        w *= (x - xi)
        xs.append(xi)
    return w, xs

def find_nodes_2(wf, a, b):
    w, xs = get_w(2)
    eq1 = integrate(expand(w * wf * qs[0]), (x, a, b))
    eq2 = integrate(expand(w * wf * qs[1]), (x, a, b))
    return solve([eq1, eq2], xs)

def find_nodes_4(wf, a, b):
    n = 4

    symbs = dict()
    poly = x**n
    for i in range(0, n):
        letter = chr(ord('A') + i)
        symbs[letter] = Symbol(letter)
        poly += (-1)**(i+1) * x**(n-i-1) * symbs[letter]

    print poly

    eqs = []
    poly *= wf
    for i in range(0, n):
        intgr = expand(poly * qs[i])
        eq = integrate(intgr, (x, a, b))
        eqs.append(eq)

    # not actual solution !
    return  solve(eqs, symbs.values())

    # w, xs = get_w(n)
    # eq1 = xs[0] + xs[1] + xs[2] + xs[3] - res[symbs['A']]
    # eq2 = xs[0]*xs[1] + xs[0]*xs[2] + xs[0]*xs[3] + xs[1]*xs[2] + xs[2]*xs[3] - res[symbs['B']]
    # eq3 = xs[0]*xs[1]*xs[2] + xs[0]*xs[1]*xs[3] + xs[0]*xs[2]*xs[3] + xs[1]*xs[2]*xs[3] - res[symbs['C']]
    # eq4 = xs[0]*xs[1]*xs[2]*xs[3] - res[symbs['D']]
    #
    # for eq in [eq1, eq2, eq3, eq4]:
    #     print eq
    #
    # return solve([eq1, eq2, eq3, eq4], *xs)

def quadrature(f, wf, nodes, a, b):
    n = len(nodes)
    w_x, xs = get_w(n)
    #print xs, nodes
    w = w_x.subs(zip(xs, nodes))
    #print w
    wdiff = w.diff(x)

    quad = 0
    AS = []
    for i in range(0, n):
        A = integrate(expand(wf * (1/wdiff.subs(x, nodes[i])) * (w/(x-nodes[i]))), (x, a, b))
        AS.append(A)
        quad += A * f.subs(x, nodes[i])

    return quad, AS




if __name__ == "__main__":

    f = cos(x)
    p = x**(-0.41)
    a = 0
    b = Rational(1, 2)

    nodes_2 = find_nodes_2(p, a, b)[0]

    nodes_4 = [0.3740045831816439, 0.02007727271691400, 0.4158706000599734, 0.13602909871867552]

    quad2, A1 = quadrature(f, p, list(nodes_2), a, b)
    quad4, A2 = quadrature(f, p, nodes_4, a, b)

    print "Integral of {}dx from 0 to 1/2".format(str(f*p))
    print ""

    print "Gauss Quadrature with 2 nodes = {:+.8f}".format(float(quad2))
    for (node, A) in zip(nodes_2, A1):
        print "xk = {:+.8f} | Ak = {:+.8f}".format(node, A)

    print ""

    print "Gauss Quadrature with 4 nodes = {:+.8f}".format(float(quad4))
    for (node, A) in zip(nodes_4, A2):
        print "xk = {:+.8f} | Ak = {:+.8f}".format(node, A)