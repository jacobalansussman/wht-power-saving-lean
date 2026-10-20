"""exact_decide.py -- w3-asym: EXACT decision of small skeletons by a Groebner basis (Singular).
Skeleton: a (2,2) core (Lemma 5 position; famI22.core22) plus the minimal landing injections the theory asks for
(x_j injected into u with a letter kappa such that D ell_e D kappa lies in output l's junk frame Phi_l, for both outputs
and both ends when found).  The unknowns are the program coefficients (merge coefficient, injection coefficients) and the
readout coefficients beta_(i,b) of every closure column, separately for each output.  Closure entries are polynomials
with integer coefficients in the program coefficients (all letters are exact group elements), so the unit equations
A_l(alpha) beta_l = e_K (l = 0, 1) form a polynomial system over Q; it has a complex solution iff its Groebner basis
is not {1} (Nullstellensatz).  Writes a Singular script per skeleton.  Usage: exact_decide.py nskel seed nsamp outdir"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from famT import *

def poly_str(polys_i, j, names):
    out = []
    for mono, vec in polys_i.items():
        c = vec[j]
        if abs(c) < 1e-12: continue
        assert abs(c.imag) < 1e-12 and abs(c.real - round(c.real)) < 1e-12
        c = int(round(c.real))
        m = "*".join(names[p] for p in mono)
        out.append(f"({c})" + (("*" + m) if m else ""))
    return "+".join(out) if out else "0"

def singular_script(P, cap=20000):
    polys = monomial_split(P)
    names = [f"a{k}" for k in range(P.np)]
    eqs = []; bvars = []
    sizes = []
    for l0 in range(P.W):
        C = Closure(P, l0, cap)
        if not C.ok: return None, None
        sizes.append((C.nr, C.nc))
        rows = {}
        tokinv = []
        for i, T in enumerate(P.terms):
            for j in range(len(T)): tokinv.append((i, j))
        for (r, c, i, j) in C.E:
            rows.setdefault(r, []).append((c, i, j))
        bn = [f"b{l0}_{c}" for c in range(C.nc)]; bvars += bn
        for r in range(C.nr):
            parts = [f"({poly_str(polys[i], j, names)})*{bn[c]}" for (c, i, j) in rows.get(r, [])]
            e = "+".join(parts) if parts else "0"
            if r == 0: e += "-1"
            eqs.append(e)
    allv = bvars + names
    s = "LIB \"poly.lib\";\n"
    s += f"ring R = 0, ({','.join(allv)}), dp;\n"
    s += "ideal I = " + ",\n  ".join(eqs) + ";\n"
    s += "option(redSB);\nideal J = std(I);\nprint(\"GB size\"); print(size(J));\nif (size(J) == 1 and J[1] == 1) { print(\"INCONSISTENT: no complex solution (not a unit)\"); } else { print(\"CONSISTENT: solutions exist; dimension\"); print(dim(J)); }\nquit;\n"
    return s, sizes

if __name__ == '__main__':
    nsk, seed, nsamp, outdir = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    for s in range(nsk):
        rec = core22(rng, 2)
        if rec is None: print("core failed"); continue
        rec2, found = guided(rng, rec, nsamp)
        # keep at most one landing letter per (output, end)
        P = build_rec(2, rec2)
        scr, sizes = singular_script(P)
        if scr is None: print(f"skel {s}: closure cap", flush=True); continue
        fn = os.path.join(outdir, f"sk{seed}_{s}.sing")
        open(fn, 'w').write(scr)
        print(f"skel {s}: injections {found}, params {P.np}, closures {sizes}, script {fn}", flush=True)
