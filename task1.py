# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 17:00:11 2014

@author: evgeniy
"""

from sympy import *

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
    z = Symbol('z')
    t = p.subs(x, 1/z) * z**p.degree()
    t = poly(simplify(t))
    if t.LC() < 0:
        t = t * (-1)
    return 1 / upper_positive(t)
    
def upper_negative(p):
    z = Symbol('z')
    t = p.subs(x, -1/z) * z**p.degree()
    t = poly(simplify(t))
    if t.LC() < 0:
        t = t * (-1)
    return -1 / upper_positive(t)
    
def lower_negative(p):
    z = Symbol('z')
    t = p.subs(x, -z) * z**p.degree()
    t = poly(simplify(t))
    if t.LC() < 0:
        t = t * (-1)
    return - upper_positive(t)

p1 = poly(x**3 + 3*x**2 - 3)

p2 = poly(2042*x**12 - 6144*x**10 + 6912*x**8 - 3524*x**6 + 840*x**4 - 72*x**2 + 1)

p3 = series(x - sin(x), x, 0, 6)

