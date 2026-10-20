"""modsolve.py (w8-b5sparse): EXACT consistency of a layered system at a random point alpha (mod p), with FLINT.

System (endsolve format): rows (non-free group elements of layer >= d), columns (y/z columns), entries (r, c, cf) with
cf = {k: v}, value = sum_k v_k a_k, a_0 = 1, a_k = alpha_k (injection / copy coefficients).  Normalisation row 0:
a random functional of the K-part of the class columns = 1 (as endsolve.assemble), or (norm='seed') the K value on the
seed columns' own K-rows only.  Consistent iff rank(A) = rank([A | e_0]).
At a uniformly random alpha in F_p^nal, consistency for generic complex alpha <=> consistency at the random point except
with probability <= deg/p (Schwartz-Zippel on the maximal minors), p = 2^61 - 1.  An INCONSISTENT verdict therefore
means: no solution for generic alpha (special alpha are not covered; endsolve's LM covers them numerically)."""
import random
from fractions import Fraction
import flint

PR = (1 << 61) - 1


def fval(v, p=PR):
    f = Fraction(v).limit_denominator(10 ** 6)
    return (f.numerator % p) * pow(f.denominator % p, p - 2, p) % p


def build(cl, alpha, rng, norm='all', seedcols=None, p=PR):
    rows, cols, ent, kl = cl
    klass, kterms = kl
    nr, nc = len(rows) + 1, len(cols)
    A = {}
    cache = {}

    def val(cf):
        key = tuple(sorted(cf.items()))
        x = cache.get(key)
        if x is None:
            x = 0
            for k, v in cf.items():
                x = (x + fval(v, p) * (1 if k == 0 else alpha[k - 1])) % p
            cache[key] = x
        return x

    for (r, c, cf) in ent:
        key = (r + 1, c)
        A[key] = (A.get(key, 0) + val(cf)) % p
    lam = {}
    use = klass if norm == 'all' else [c for c in klass if seedcols is None or c in seedcols]
    seedrows = None
    if norm == 'seed' and seedcols is not None:
        seedrows = set(rk for c in seedcols for (rk, cf) in kterms.get(c, []))
    for c in klass:
        for (rk, cf) in kterms.get(c, []):
            if seedrows is not None and rk not in seedrows:
                continue
            if rk not in lam:
                lam[rk] = rng.randrange(1, p)
            key = (0, c)
            A[key] = (A.get(key, 0) + lam[rk] * val(cf)) % p
    return nr, nc, A


def consistent(cl, nal, seed=1, norm='all', seedcols=None, p=PR):
    rng = random.Random(seed)
    alpha = [rng.randrange(1, p) for _ in range(nal)]
    nr, nc, A = build(cl, alpha, rng, norm, seedcols, p)
    if nc == 0:
        return False, 0, 0, nr, nc
    Ma = flint.nmod_mat(nr, nc, p)
    Mt = flint.nmod_mat(nr, nc + 1, p)
    for (r, c), v in A.items():
        if v:
            Ma[r, c] = v; Mt[r, c] = v
    Mt[0, nc] = 1
    ra = Ma.rank(); rt = Mt.rank()
    return (ra == rt), ra, rt, nr, nc
