# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 17:00:11 2014

@author: evgeniy
"""

from sympy import *
from sympy.abc import x, z

def upper_positive(p):
    coeffs = p.all_coeffs()
    i = 0
    m = 0
    a = 0
    a_0 = coeffs[0]
    for coeff in coeffs:
        if coeff < 0:
            m = i
        if abs(coeff) > a:
            a = abs(coeff)
        i += 1
    return 1 + (a / a_0)**(1/m)
    
def lower_positive(p):
    t = p.subs(x, 1/z) * z**p.degree()
    t = poly(simplify(t))
    if t.LC() < 0:
        t = t * (-1)
    return 1 / upper_positive(t)
    
def upper_negative(p):
    t = p.subs(x, -1/z) * z**p.degree()
    t = poly(simplify(t))
    if t.LC() < 0:
        t = t * (-1)
    return -1 / upper_positive(t)
    
def lower_negative(p):
    t = p.subs(x, -z) * z**p.degree()
    t = poly(simplify(t))
    if t.LC() < 0:
        t = t * (-1)
    return - upper_positive(t)

def f_range(a, b, d):
    while a < b:
        yield a
        a += d

def find_intervals(a, b):
    if a > b:
        tmp = a
        a = b
        b = tmp
    intervals = list()
    j = a    
    for i in f_range(a, b, 0.1):
        if (i * j) < 0:
            intervals.append((i, j))
        j = i
    return intervals
    
def find_init_approx(p, interval):
    for x_0 in f_range(interval[0], interval[1], 0.01):
        if abs(p.subs(x, x_0)) < 0.1:
            return x_0

def newtone(x_0, p):
    k = 1
    x_k = x_0
    while k < 10:
        x_k = x_k - p.subs(x, x_k) / diff(p, x).subs(x, x_k)
        print "{}    |    {}    |    {}".format(k, x_k, p.subs(x, x_k))
   
p1 = poly(x**3 + 3*x**2 - 3)

p2 = poly(2042*x**12 - 6144*x**10 + 6912*x**8 - 3524*x**6 + 840*x**4 - 72*x**2 + 1)

p3 = series(x - sin(x), x, 0, 6)
     
def solve(p):
    k_1 = upper_positive(p)
    k_2 = lower_positive(p)
    k_3 = lower_negative(p)
    k_4 = upper_negative(p)
    intrvs = find_intervals(k_1, k_2) + find_intervals(k_3, k_4)
    for intrv in intrvs:
        x_0 = find_init_approx(p, intrv)
        print "-------------------------------------------------------------------------------"
        newtone(x_0, p)

if __name__ == "__main__":
    solve(p1)

