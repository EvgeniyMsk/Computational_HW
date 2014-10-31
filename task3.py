from sys import stdout
from sympy import *

def f_range(a, b, d):
    while a <= b:
        yield a
        a += d

if __name__ == "__main__":

    n = 10

    h = 0.1

    xs   = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    f_xs = [1.614419,
            1.656832,
            1.694888,
            1.728606,
            1.758030,
            1.783225,
            1.804279,
            1.821299,
            1.834414,
            1.843768]

    #print zip(xs, f_xs)

    m = 5

    diffs = []
    diffs.append(f_xs)

    # calculate differences
    for i in range(1, m):
        i_diff      = []
        prev_diff   = diffs[i - 1]
        for j in range(1, len(prev_diff)):
            i_diff.append(prev_diff[j] - prev_diff[j - 1])
        diffs.append(i_diff)

    lens = list()
    for i in range(0, len(diffs)):
        max = 0
        for j in range(0, len(diffs[i])):
            if max < len(str(diffs[i][j])):
                max = len(str(diffs[i][j]))
        lens.append(max + 3)

    # print differences
    for i in range(0, len(diffs[0])):
        stdout.write(str(xs[i]) + '  ')
        for j in range(0, len(diffs)):
            if i < len(diffs[j]):
                stdout.write(str(diffs[j][i]) + ' '*(lens[j] - len(str(diffs[j][i]))))
        stdout.write('\n')

    # interpolation polynomial
    t  = Symbol('t')
    y0 = Symbol('y0')
    y1 = Symbol('y1')
    y2 = Symbol('y2')
    y3 = Symbol('y3')
    y4 = Symbol('y4')

    intrp1 = Poly(y0 + y1 * t + (y2 / 2) * t * (t - 1) + (y3 / 6) * t * (t - 1) * (t - 2)
                     + (y4 / 24) * t * (t - 1) * (t - 2) * (t - 3))

    intrp2 = Poly(y0 + y1 * t + (y2 / 2) * t * (t + 1) + (y3 / 6) * t * (t + 1) * (t + 2)
                     + (y4 / 24) * t * (t + 1) * (t + 2) * (t + 3))

    intrp3 = Poly(y0 + y1 * t + (y2 / 2) * t * (t - 1) + (y3 / 6) * t * (t - 1) * (t + 1)
                     + (y4 / 24) * t * (t - 1) * (t + 1) * (t - 2))


    x1 = 0.171494
    x2 = 0.764488
    x3 = 0.408541

    y = 1.817978

    sub1 = [(y0, diffs[0][1]),
            (y1, diffs[1][1]),
            (y2, diffs[2][1]),
            (y3, diffs[3][1]),
            (y4, diffs[4][1]),
            (t,  (x1 - xs[1]) / h)]

    sub2 = [(y0, diffs[0][n - 2]),
            (y1, diffs[1][n - 3]),
            (y2, diffs[2][n - 4]),
            (y3, diffs[3][n - 5]),
            (y4, diffs[4][n - 6]),
            (t,  (x2 - xs[8]) / h)]


    p1 = intrp1.subs(sub1)
    p2 = intrp2.subs(sub2)


    print " "

    print "       x_1      |       p(x_1)   "
    print "---------------------------------------------------------------------------"
    print "  {:+.9f}  |  {:+.9f}".format(float(x1), float(p1))
    print "---------------------------------------------------------------------------"

    print "       x_2      |       p(x_2)   "
    print "---------------------------------------------------------------------------"
    print "  {:+.9f}  |  {:+.9f}".format(float(x2), float(p2))
    print "---------------------------------------------------------------------------"


    ind = 4
    sub1 = [(y0, diffs[0][ind]),
            (y1, diffs[1][ind]),
            (y2, diffs[2][ind]),
            (y3, diffs[3][ind]),
            (y4, diffs[4][ind]),
            (t,  (x3 - xs[ind]) / h)]
    p3 = intrp1.subs(sub1)

    print
    print " Table begin "
    print
    print "       x_3      |       p(x_3)   "
    print "---------------------------------------------------------------------------"
    print "  {:+.9f}  |  {:+.9f}".format(float(x3), float(p3))
    print "---------------------------------------------------------------------------"


    sub3 = [(y0, diffs[0][ind]),
            (y1, diffs[1][ind]),
            (y2, diffs[2][ind - 1]),
            (y3, diffs[3][ind - 1]),
            (y4, diffs[4][ind - 2]),
            (t,  (x3 - xs[ind]) / h)]

    p3 = intrp3.subs(sub3)

    print
    print " Table middle "
    print
    print "       x_3      |       p(x_3)   "
    print "---------------------------------------------------------------------------"
    print "  {:+.9f}  |  {:+.9f}".format(float(x3), float(p3))
    print "---------------------------------------------------------------------------"

    sub3 = [(y0, diffs[0][n - ind - 1]),
            (y1, diffs[1][n - ind - 2]),
            (y2, diffs[2][n - ind - 3]),
            (y3, diffs[3][n - ind - 4]),
            (y4, diffs[4][n - ind - 5]),
            (t,  (x3 - xs[n - ind - 1]) / h)]

    p3 = intrp2.subs(sub3)

    print
    print " Table end "
    print
    print "       x_3      |       p(x_3)   "
    print "---------------------------------------------------------------------------"
    print "  {:+.9f}  |  {:+.9f}".format(float(x3), float(p3))
    print "---------------------------------------------------------------------------"

    ind = 6
    sub4 = [(y0, diffs[0][ind]),
            (y1, diffs[1][ind]),
            (y2, diffs[2][ind - 1]),
            (y3, diffs[3][ind - 1]),
            (y4, diffs[4][ind - 2])]
    intrpr4 = intrp3
    p_t = ((intrpr4 - y - t*y1) / y1).subs(sub4)



    print
    print " Inverse interpolation "
    print

    ind = 6
    sub4 = [(y0, diffs[0][ind]),
            (y1, diffs[1][ind]),
            (y2, diffs[2][ind - 1]),
            (y3, diffs[3][ind - 1]),
            (y4, diffs[4][ind - 2])]
    intrpr4 = intrp3
    p_t = ((intrpr4 - y - t*y1) / y1).subs(sub4)

    tk   = 0
    tk_1 = 0
    eps = 1e-5
    for k in range(0, 20):
        tk_1 = p_t.subs(t, tk)
        print " k  |       tk       |       p(tk)    "
        print "---------------------------------------------------------------------------"
        print " {}   |  {:+.9f}  |  {:+.9f}".format(k, float(tk), float(tk_1))
        print "---------------------------------------------------------------------------"
        if abs(tk_1 - tk) < eps:
            break
        tk = tk_1

    tk = tk_1
    x = h * tk + xs[ind]

    p_x = intrpr4.subs(sub4)
    p_x = p_x.subs(t, tk)

    print ""
    print "       y*       |       x*       |       p(x*)   "
    print "---------------------------------------------------------------------------"
    print "  {:+.9f}  |  {:+.9f}  |  {:+.9f}".format(float(y), float(x), float(p_x))
    print "---------------------------------------------------------------------------"

