# vc_run.py -- v-cert-07 driver.  usage:
#   vc_run.py hook                 -> unitarity, invariance, hook psi/phi (float, complex d) + psi mod 2 primes
#   vc_run.py cuts MODE b1,b2,...  -> psi_1, psi_2, psi_3 (float, MODE f = complex d, r = real rotation) of orbit-union cuts
#                                     (bits in cert-07's orbit numbering, reproduced independently here)
#   vc_run.py cover                -> orbit sizes, mirror map, number of canonical cuts
import sys, time, json
import numpy as np
from vclib import *
t0 = time.time()
P1, P2 = 998244353, 1000000009
def T(): return f"({time.time()-t0:.0f}s)"
orb, lab, gens = orbits_gl42()
K = len(orb)
print("orbit sizes", [len(o) for o in orb], "K", K, T(), flush=True)
hook_orbits = sorted(set(int(x) for x in lab[HOOK]))
print("hook orbits", hook_orbits, "dim", int(HOOK.sum()), "pure union:", all(HOOK[o].all() or (~HOOK[o]).all() for o in orb), flush=True)
def inS_of(b): return np.array([(b >> int(lab[k])) & 1 == 1 for k in range(NT)])
mode = sys.argv[1]
if mode == 'cover':
    comp = np.array([IDX[tuple(sorted(set(range(16)) - set(s)))] for s in SUBS])
    mir = [int(lab[comp[o[0]]]) for o in orb]
    assert all(len(set(lab[comp[o]])) == 1 for o in orb)
    full = (1 << K) - 1
    def canon(b):
        c = []
        for x in (b, sum(1 << mir[k] for k in range(K) if b >> k & 1)):
            if x >> (K - 1) & 1: x = full ^ x
            c.append(x)
        return min(c)
    can = sorted({canon(b) for b in range(1, 1 << (K - 1))})
    print("mirror", mir, "canonical cuts", len(can), flush=True)
    json.dump(can, open("canon.json", "w"))
    sys.exit()
if mode == 'hook':
    for j in (1, 2, 3):
        B = blocks(j, 'f')
        err = max(np.abs(G @ G.conj().T - np.eye(len(r))).max() for r, G in B.values())
        print(f"j={j}: {len(B)} classes, max {max(len(r) for r, G in B.values())}, unitarity err {err:.1e}", T(), flush=True)
        ps, ph, lo, hi = psi_float(B, HOOK)
        print(f"  HOOK float: psi_{j} = {ps} (x16 = {16*ps}), phi_{j} = {ph:.6f} (x16 = {16*ph:.3f}), s_{j} = {ph/ps:.6f}, kept >= {lo:.3e}, dropped <= {hi:.1e}", T(), flush=True)
        # invariance under a generic diagonal Clifford letter: rho(D)[T',T] = D[sT',sT] det D[T',T]/D00 is diagonal
    # diagonal letters: check rho of D = diag(i^q(x)) with q = random quadratic form is diagonal & the hook is a union of orbits
    rng = np.random.default_rng(7)
    q = rng.integers(0, 4, 16); Dd = np.diag(1j ** q)
    vals = [Dd[SIG[k], SIG[k]] * np.prod([Dd[t, t] for t in SUBS[k]]) / Dd[0, 0] for k in range(5)]
    print("diagonal letter acts on f_T by scalars (structural), sample", np.round(vals, 3), flush=True)
    for p in (P1, P2):
        for j in (1, 2, 3):
            B = blocks(j, 'p', p)
            print(f"  HOOK mod {p}: psi_{j} = {psi_modp(B, HOOK, p)}", T(), flush=True)
    sys.exit()
if mode == 'cuts':
    md = sys.argv[2]; bl = [int(x) for x in sys.argv[3].split(',')]
    BB = {j: blocks(j, md) for j in (1, 2, 3)}
    print("blocks built", T(), flush=True)
    for b in [1433] + bl:
        inS = inS_of(b)
        out = dict(bits=b, orbits=[k for k in range(K) if b >> k & 1], dimN=int(inS.sum()), ub3=ub3(BB[3], inS))
        for j in (1, 2, 3):
            ps, ph, lo, hi = psi_float(BB[j], inS)
            out[f"psi{j}"] = ps; out[f"phi{j}"] = round(ph, 6); out[f"gap{j}"] = (lo, hi)
        out["ratio3"] = out["psi3"] / out["psi1"]; out["spec3"] = out["phi3"] / out["phi1"]
        out["margin_vs_hook"] = out["psi3"] - 2135 / 858 * out["psi1"]; out["margin_vs_2.5"] = out["psi3"] - 2.5 * out["psi1"]
        print(json.dumps(out), T(), flush=True)
