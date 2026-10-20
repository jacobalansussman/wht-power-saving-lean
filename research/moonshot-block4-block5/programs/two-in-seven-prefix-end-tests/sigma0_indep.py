"""sigma0_indep.py -- o3-finish2: a second implementation of Criterion 0 (own rebuild with plain numpy matrices, own
K-frame test, bm's random-matching decomposition fulltestG.decomp instead of kfG.decomp1, own frame canonicalisation
from frametest.frame_data for the class key).  Compares dim0 with out_sigma0_all.jsonl.
usage: python3 sigma0_indep.py PREFIXES.jsonl SIGMA0_OUT.jsonl"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
PYLIB = ''
sys.path.insert(0, PYLIB)
import numpy as np
from fractions import Fraction
sys.path.insert(0, HERE)
from frametest import frame_data, toint
sys.path.insert(1, HERE)   # fulltestG.py sits in this folder
import fulltestG as FG

A, B = (1 + 1j) / 2, (1 - 1j) / 2
n = 16
I = np.eye(n, dtype=complex)
def X(z):
    M = np.zeros((n, n), dtype=complex)
    for x in range(n):
        M[x ^ z, x] = 1
    return M
D = {z: A * I + B * X(z) for z in (1, 2, 4, 8)}
K = I.copy()
for j in range(4):
    K = K @ D[1 << j]
Z4 = [1, 1j, -1, -1j]

def mono(s, q):
    M = np.zeros((n, n), dtype=complex)
    for y in range(n):
        M[s[y], y] = Z4[q[y] % 4]
    return M

def key(M):
    re, im = toint(M)
    return re.tobytes() + im.tobytes()

def lmul(M, gen):
    out = {}
    for (l, k), (G, c) in gen.items():
        H = M @ G
        kk = (l, key(H))
        if kk in out:
            out[kk] = (out[kk][0], out[kk][1] + c)
        else:
            out[kk] = (H, c)
    return {k: v for k, v in out.items() if v[1] != 0}

def add(a, b, cb):
    out = dict(a)
    for k, (G, c) in b.items():
        if k in out:
            out[k] = (out[k][0], out[k][1] + cb * c)
        else:
            out[k] = (G, cb * c)
    return {k: v for k, v in out.items() if v[1] != 0}

def build(prog):
    gens = [{(0, key(I)): (I, Fraction(1))}, {(1, key(I)): (I, Fraction(1))}]
    for mv in prog:
        b = mono(*mv["b"])
        c = {}
        for (nu, de, a, j) in mv["cs"]:
            t = gens[j] if a is None else lmul(mono(*a), gens[j])
            c = add(c, t, Fraction(nu, de))
        gens.append(lmul(D[1] @ b, c))
    return gens

def is_mono(M):
    nz = np.abs(M) > 1e-9
    if not ((nz.sum(axis=0) == 1).all() and (nz.sum(axis=1) == 1).all()):
        return False
    return all(min(abs(v - z) for z in Z4) < 1e-9 for v in M[nz])

Kinv = np.conj(K).T
def in_kframe(G):
    return is_mono(G @ Kinv)

class _G:          # the minimal Game interface fulltestG.decomp needs
    import bm as _bm
    def __init__(self):
        self.g = self._bm.Game(4)
        self.D1 = self.g.D1
    def __getattr__(self, a):
        return getattr(self.g, a)

def sig(gen):
    s = [Fraction(0), Fraction(0)]
    for (l, k), (G, c) in gen.items():
        if in_kframe(G):
            s[l] += c
    return s

def rank2(vs):
    vs = [v for v in vs if v != [0, 0]]
    if not vs:
        return 0
    a = vs[0]
    return 2 if any(a[0] * b[1] - a[1] * b[0] != 0 for b in vs[1:]) else 1

if __name__ == "__main__":
    GG = _G()
    ref = {}
    for l in open(sys.argv[2]):
        r = json.loads(l); ref[(r["run"], r["it"])] = r["dim0"]
    agree = disagree = fails = 0
    for l in open(sys.argv[1]):
        rec = json.loads(l)
        k = (rec["run"].replace(".jsonl", ""), rec["it"])
        gens = build(rec["prog"])
        vs = []
        for t in range(2, len(gens)):
            if any(in_kframe(G) for (G, c) in gens[t].values()):
                vs.append(sig(gens[t]))
        seen = set()
        for s in range(2, len(gens)):
            for (l2, kk), (G, c) in gens[s].items():
                h = GG.mul(K, np.conj(G).T)
                if not (np.abs(h) > 1e-9).sum(axis=1).max() == 2:
                    continue
                d = FG.decomp(GG, h)
                if d is None:
                    if (np.abs(h) > 1e-9).sum(axis=1).min() == 2:
                        fails += 1
                    continue
                b3, b4 = d
                ck = (s, frame_data(D[1] @ b4)[0])
                if ck in seen:
                    continue
                seen.add(ck)
                vs.append(sig(lmul(D[1] @ b4, gens[s])))
        d0 = rank2(vs)
        if d0 == ref.get(k):
            agree += 1
        else:
            disagree += 1
            print("DISAGREE", k, d0, ref.get(k), flush=True)
    print("agree", agree, "disagree", disagree, "decomposition failures on support-2 elements", fails)
