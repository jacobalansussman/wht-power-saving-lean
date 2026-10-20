"""exact.py -- w5-twoshared: exact decision of fixed-letter skeletons of the two-shared shape, with EVERY coefficient
free (meta2 FOLLOWUP notice: fixing a first summand to 1 drops the zero pattern).  Closure = asym.Closure (bm.Game.closure,
exact over all f); polynomial system = closure equations of every output with the readout coefficient of every closure
column as unknowns; Groebner basis {1} <=> no complex coefficients make the skeleton a unit (Nullstellensatz)."""
import sys, os, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
import numpy as np
from asym import G, Prog, Closure, monomial_split

class ProgF(Prog):
    """every summand of every move carries its own free coefficient"""
    def add_move(self, summ, D=None):
        D = G.D1 if D is None else D
        terms = {}; order = []; maps = []; rec = []
        for k, (src, b) in enumerate(summ):
            pidx = self.np + k
            si, ri = [], []
            for j, (l, g) in enumerate(self.terms[src]):
                g3 = G.intern(G.mul(D, b, G.MATS[g]))
                key = (l, g3)
                if key not in terms:
                    terms[key] = len(order); order.append(key)
                si.append(j); ri.append(terms[key])
            maps.append((pidx, src, np.array(si), np.array(ri)))
            rec.append((src, b, pidx))
        self.np += len(summ)
        self.terms.append(order); self.maps.append(maps); self.moves.append(rec)
        return len(self.terms) - 1

def build(W, moves):
    """moves: list of summand lists [(src, b)], src = generator index (inputs 0..W-1)"""
    P = ProgF(W)
    for summ in moves:
        P.add_move(summ)
    return P

def poly_str(polys_i, j, names):
    out = []
    for mono, vec in polys_i.items():
        c = vec[j]
        if abs(c) < 1e-12: continue
        assert abs(c.imag) < 1e-12 and abs(c.real - round(c.real)) < 1e-12, c
        c = int(round(c.real))
        m = "*".join(names[p] for p in mono)
        out.append(f"({c})" + (("*" + m) if m else ""))
    return "+".join(out) if out else "0"

def singular_script(P, cap=20000, outputs=None):
    polys = monomial_split(P)
    names = [f"a{k}" for k in range(P.np)]
    eqs = []; bvars = []; sizes = []
    outputs = range(P.W) if outputs is None else outputs
    for l0 in outputs:
        C = Closure(P, l0, cap)
        if not C.ok: return None, None
        sizes.append((C.nr, C.nc))
        rows = {}
        for (r, c, i, j) in C.E:
            rows.setdefault(r, []).append((c, i, j))
        bn = [f"b{l0}_{c}" for c in range(C.nc)]; bvars += bn
        for r in range(C.nr):
            parts = [f"({poly_str(polys[i], j, names)})*{bn[c]}" for (c, i, j) in rows.get(r, [])]
            e = "+".join(parts) if parts else "0"
            if r == 0: e += "-1"
            eqs.append(e)
    allv = bvars + names
    s = f"ring R = 0, ({','.join(allv)}), dp;\n"
    s += "ideal I = " + ",\n  ".join(eqs) + ";\n"
    s += "option(redSB);\nideal J = std(I);\nprint(\"GB size\"); print(size(J));\n"
    s += "if (size(J) == 1 and J[1] == 1) { print(\"VERDICT INCONSISTENT\"); } else { print(\"VERDICT CONSISTENT dim\"); print(dim(J)); }\nquit;\n"
    return s, sizes

def decide(scr, fn, timeout=120):
    open(fn, 'w').write(scr)
    try:
        out = subprocess.run(['Singular', '-q', fn], capture_output=True, text=True, timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        return 'TIMEOUT', ''
    if 'VERDICT INCONSISTENT' in out: return 'INCONSISTENT', out
    if 'VERDICT CONSISTENT' in out: return 'CONSISTENT', out
    return 'ERR', out
