"""twist_check.py (v-ex2-limits, 23:20): in the ENGINE game the letter group contains the scalar i*I, so the f-limit
along one ultrafilter only sees one central character (residue of f mod 4).  Pure python, exact Gaussian rationals."""
from fractions import Fraction as Fr
def cm(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def ca(a, b): return (a[0]+b[0], a[1]+b[1])
def mm(A, B):
    n = len(A); return [[(lambda s: s)(sum_c([cm(A[i][k], B[k][j]) for k in range(n)])) for j in range(n)] for i in range(n)]
def sum_c(l):
    s = (Fr(0), Fr(0))
    for x in l: s = ca(s, x)
    return s
a = (Fr(1, 2), Fr(1, 2)); b = (Fr(1, 2), Fr(-1, 2)); z = (Fr(0), Fr(0)); one = (Fr(1), Fr(0)); i = (Fr(0), Fr(1))
DIR = [[a, b], [b, a]]          # dir(e_1) = a I + b X   (paid letter, block 1)
S = [[one, z], [z, i]]          # S = diag(1, i)         (Clifford Z4 phase, free in E)
W = mm(mm(S, DIR), S)           # S.dir.S = omega * H
W2 = mm(W, W)
print("S.dir.S =", W)
print("(S.dir.S)^2 =", W2, " -> equals i*I:", W2 == [[i, z], [z, i]])
# family of the word (S.D.S)^2 at f columns is (i I)^(x)f = i^f I, so x = 1 - (S.D.S)^2 has x_f = (1 - i^f) I:
for f in range(1, 9):
    s = one
    for _ in range(f): s = cm(s, i)
    r = ca(one, (-s[0], -s[1])); rho = r[0]**2 + r[1]**2
    print("f=%d: rho_x(f) = |1 - i^f|^2 = %s" % (f, rho))
