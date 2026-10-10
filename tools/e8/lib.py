#!/usr/bin/env python3
"""lib.py -- E8 family helpers for the search programs.  Standard library only."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import e8
from fractions import Fraction as Fr
P = e8.PORTS; V = 120; IDX = e8.INDEX; G = e8.GRAM; ONE = 511
def par(x): return bin(x).count("1") & 1
def transvect(w):
    """the reflection in the root of port w, as a map on labels: u -> u + (u.w')w', w' = w + all-ones"""
    wp = w ^ ONE
    return lambda u: u ^ wp if par(u & wp) else u
def perm_map(pi):
    """point permutation pi (list: i -> pi[i]) acting on labels"""
    def f(u):
        v = 0
        for i in range(9):
            if u >> i & 1: v |= 1 << pi[i]
        return v
    return f
def compose(*fs):
    def f(u):
        for g in reversed(fs): u = g(u)
        return u
    return f
def as_perm(f):
    """the permutation of port indices induced by a label map (None if it does not permute the ports)"""
    try:
        p = [IDX[f(u)] for u in P]
    except KeyError:
        return None
    return p if len(set(p)) == V else None
def omega():
    """an order-3 map without fixed ports: point permutation (012)(345)(678) times the rotation of the A2 spanned
    by the triples 012, 345, 678 (product of two reflections)"""
    pi = perm_map([1, 2, 0, 4, 5, 3, 7, 8, 6])
    t1 = transvect(0b000000111); t2 = transvect(0b000111000)
    return compose(pi, t1, t2)
def need(S):
    """what y_S holds too much after the scatter with the matrix B: {T: B[S][T]} for T != S"""
    return {T: Fr(G[S][T], 2) for T in range(V) if T != S and G[S][T]}
