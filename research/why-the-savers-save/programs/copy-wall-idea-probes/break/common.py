"""common.py (g-break) -- E8 lines, F_2 labels u_S (9 bits, u_S.u_T = B0 = 1+G mod 2), helpers.  Imports u-floor/copyfloor."""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-floor-checker'))
import copyfloor as cf

L, G = cf.e8_lines()
B0 = [[(1 + G[i][j]) % 2 for j in range(120)] for i in range(120)]
U = cf.f2_labels(B0)
UI = [sum(U[i][c] << c for c in range(9)) for i in range(120)]   # labels as 9-bit ints

def f2span(xs):
    basis = []
    for x in xs:
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x); basis.sort(reverse=True)
    return basis

def f2dim(xs):
    return len(f2span(xs))

def dot9(a, b):
    return bin(a & b).count('1') & 1
