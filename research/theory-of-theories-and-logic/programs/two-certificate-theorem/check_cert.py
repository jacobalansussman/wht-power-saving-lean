"""check_cert.py -- ml-fooling: STANDALONE checker for commuting-matrix existence certificates (NOTE.md s2.5).
Imports nothing from this folder; uses only Python's fractions and json.
A certificate file (JSON) holds:
  "unknowns": list of names (n);
  "equations": list of polynomials, each a list of [exponent list (length n), "p/q" coefficient];
  "matrices": list of n matrices, each r x r, entries "p/q";
  "provenance": free text (how the system was sliced; NOT used by the check).
The check: (a) the n matrices commute pairwise; (b) every equation evaluated at the matrices is the zero matrix.
If both hold, the equations have a common complex zero (a common eigenvector of the matrices gives one). For a sliced
system, a zero of the slice extends to a zero of the original system through the substitution in "provenance".
usage: check_cert.py CERT.json [...]"""
import sys, json
from fractions import Fraction

def check(path):
    C = json.load(open(path))
    n = len(C['unknowns']); Ms = [[[Fraction(x) for x in row] for row in M] for M in C['matrices']]
    assert len(Ms) == n, 'need one matrix per unknown'
    r = len(Ms[0]); assert r >= 1 and all(len(M) == r and all(len(row) == r for row in M) for M in Ms), 'shapes'
    def mm(A, B): return [[sum(A[i][t] * B[t][j] for t in range(r)) for j in range(r)] for i in range(r)]
    I = [[Fraction(int(i == j)) for j in range(r)] for i in range(r)]
    for i in range(n):
        for j in range(i + 1, n):
            if mm(Ms[i], Ms[j]) != mm(Ms[j], Ms[i]): return False, f'matrices {i} and {j} do not commute'
    pw = {}
    def mpow(i, k):
        if k == 0: return I
        if (i, k) not in pw: pw[(i, k)] = mm(mpow(i, k - 1), Ms[i])
        return pw[(i, k)]
    for q, f in enumerate(C['equations']):
        S = [[Fraction(0)] * r for _ in range(r)]
        for e, c in f:
            assert len(e) == n
            A = I
            for i, k in enumerate(e):
                if k: A = mm(A, mpow(i, k))
            c = Fraction(c)
            S = [[S[a][b] + c * A[a][b] for b in range(r)] for a in range(r)]
        if any(x != 0 for row in S for x in row): return False, f'equation {q} does not vanish'
    msg = f'{n} unknowns, {len(C["equations"])} equations, {r} x {r} matrices: commute and all equations vanish'
    if 'original' in C:
        # the slice must be the original system with each eliminated unknown replaced by an affine form in the kept ones
        on = C['original']['unknowns']; keep = C['unknowns']; sub = C['substitution']
        assert set(keep) | set(sub) == set(on) and not (set(keep) & set(sub)), 'unknowns do not partition'
        kpos = {v: t for t, v in enumerate(keep)}
        forms = {v: {tuple(e): Fraction(c) for e, c in f} for v, f in sub.items()}
        assert all(sum(e) <= 1 for f in forms.values() for e in f), 'substitution is not affine'
        def pmul(a, b):
            out = {}
            for ea, ca in a.items():
                for eb, cb in b.items():
                    k = tuple(x + y for x, y in zip(ea, eb)); out[k] = out.get(k, 0) + ca * cb
            return {k: v for k, v in out.items() if v}
        z = tuple([0] * len(keep))
        def canon(eqs): return sorted(sorted((tuple(e), Fraction(c)) for e, c in f.items()) for f in eqs)
        got = []
        for f in C['original']['equations']:
            acc = {}
            for e, c in f:
                term = {z: Fraction(c)}
                for v, k in zip(on, e):
                    for _ in range(k):
                        if v in kpos:
                            u = [0] * len(keep); u[kpos[v]] = 1; term = pmul(term, {tuple(u): Fraction(1)})
                        else: term = pmul(term, forms[v])
                for kk, vv in term.items(): acc[kk] = acc.get(kk, 0) + vv
            acc = {k: v for k, v in acc.items() if v}
            if acc: got.append(acc)
        want = [{tuple(e): Fraction(c) for e, c in f} for f in C['equations']]
        if canon(got) != canon(want): return False, msg + '; BUT the slice is not the original system under the substitution'
        msg += f'; slice = original system ({len(on)} unknowns, {len(C["original"]["equations"])} equations) under an affine substitution of {len(sub)} unknowns, so the original has a complex zero'
    return True, msg

if __name__ == '__main__':
    for p in sys.argv[1:]:
        ok, msg = check(p)
        print(('VALID   ' if ok else 'INVALID ') + p + ': ' + msg, flush=True)
