"""annih2.py -- ms-xw3.  Cross-talk test, version 2: identities required at TWO f at once (rows of f1 and f2 stacked,
one coefficient vector), then each found identity re-tested at a third f (f3) by direct evaluation.

Question: are there bounded mixes  sum_{g,h} c_{gh} g m h (y)  (y = M1 G M0 v, a private 2-layer chain; m one paid
move) that equal a combination of <= 2-layer private terms  sum_b (a_b b(y) + a'_b b(M0 v) + a''_b b(v)),  with the
combination NOT zero itself?  Such an identity is what a shared top move needs (NOTE section 3).
usage: python3 -I annih2.py <menuR> <menuC> <menuB> <chain> <z_m[,C_m]> <f1,f2> <f3> [S]
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw3lib import Game, menu as MENU


def chain_vecs(G, chain):
    q = chain.split(";")
    if chain == "id":
        v = G.inp(0); return v, v, v
    if len(q) == 1:                          # one private move: y = M0 v; returns (v, v, y)
        z0, C0 = (q[0].split(",") + ["1"])[:2]
        v = G.inp(0)
        return v, v, G.move(int(z0), C0, v)
    z0, C0 = (q[0].split(",") + ["1"])[:2]; Gw = q[1]; z1, C1 = (q[2].split(",") + ["1"])[:2]
    v = G.inp(0)
    y0 = G.move(int(z0), C0, v)
    y1 = G.move(int(z1), C1, G.W_(Gw, y0))
    return v, y0, y1


def fl(G, X):
    return np.concatenate([X[f].reshape(-1) for f in G.FS])


def build(G, MR, MC, MB, chain, zm, Cm):
    v, y0, y1 = chain_vecs(G, chain)
    A3 = []
    for h in MC:
        u = G.move(zm, Cm, G.W_(h, y1))
        for g in MR: A3.append(fl(G, G.W_(g, u)))
    A2 = [fl(G, G.W_(b, s)) for s in (y1, y0, v) for b in MB]
    return np.stack(A3, 1), np.stack(A2, 1)


def deg_nontrivial(A, Ap, triv):
    """dimension of span(non-trivial columns) cap L2, from the full and projected column sets"""
    import numpy as _np
    keep = _np.setdiff1d(_np.arange(A.shape[1]), triv)
    if keep.size == 0: return 0
    s1 = _np.linalg.svd(A[:, keep], compute_uv=False); r1 = int((s1 > s1[0] * 1e-9).sum())
    s2 = _np.linalg.svd(Ap[:, keep], compute_uv=False); r2 = int((s2 > s2[0] * 1e-9).sum())
    return r1 - r2


def main():
    mR, mC, mB, chain, zmC, fs, f3 = sys.argv[1:8]
    S = int(sys.argv[8]) if len(sys.argv) > 8 else 1
    zm, Cm = (zmC.split(",") + ["1"])[:2]; zm = int(zm)
    MR, MC, MB = MENU(mR), MENU(mC), MENU(mB)
    FS = [int(a) for a in fs.split(",")]
    t0 = time.time()
    G = Game(FS, W=1, S=S, seed=5)
    A3, A2 = build(G, MR, MC, MB, chain, zm, Cm)
    # orthonormal basis of span A2, project it out of A3: degenerate combos = null space of (I - P2) A3 minus null of A3
    U2, s2, _ = np.linalg.svd(A2, full_matrices=False)
    Q2 = U2[:, s2 > s2[0] * 1e-9]
    A3p = A3 - Q2 @ (Q2.conj().T @ A3)
    sA = np.linalg.svd(A3, compute_uv=False); rA = int((sA > sA[0] * 1e-9).sum())
    U, sP, Vh = np.linalg.svd(A3p, full_matrices=False); rP = int((sP > sP[0] * 1e-9).sum())
    inter = rA - rP
    # columns that are individually <= 2-layer (the trivial collapse of a repeated direction): their span is "trivial"
    cn = np.linalg.norm(A3p, axis=0) / np.maximum(np.linalg.norm(A3, axis=0), 1e-300)
    triv = np.flatnonzero(cn < 1e-8)
    rT = 0
    if triv.size:
        sT = np.linalg.svd(A3[:, triv], compute_uv=False); rT = int((sT > sT[0] * 1e-9).sum())
    excess = inter - rT
    print(f"# annih2 R={mR}({len(MR)}) C={mC}({len(MC)}) B={mB}({len(MB)}) chain={chain} m=({zm},{Cm}) FS={FS} rows={A3.shape[0]} "
          f"rank A3 {rA}, rank A3 mod L2 {rP}: degenerate dimension {inter}, trivial {rT}, EXCESS {excess}, NONTRIVIAL-COLLAPSE {deg_nontrivial(A3, A3p, triv)}  [sv near cut: {sP[max(rP-2,0):rP+2]}] ({time.time()-t0:.0f}s)",
          flush=True)
    if excess <= 0 and os.environ.get('SHOWALL') != '1': return
    # coefficient vectors x with A3p x = 0 but A3 x != 0
    N = Vh[rP:].conj().T
    W = A3 @ N
    u, sw, vh = np.linalg.svd(W, full_matrices=False)
    kgood = int((sw > sw[0] * 1e-8).sum()) if sw.size else 0
    Xg = N @ vh[:kgood].conj().T
    G3 = Game([int(f3)], W=1, S=1, seed=6)
    v3, y03, y13 = chain_vecs(G3, chain)
    B2 = np.stack([fl(G3, G3.W_(b, s)) for s in (y13, y03, v3) for b in MB], 1)
    for k in range(min(kgood, 10)):
        x = Xg[:, k]
        lhs = np.zeros(B2.shape[0], complex)
        for hi, h in enumerate(MC):
            cs = x[hi * len(MR):(hi + 1) * len(MR)]
            if np.abs(cs).max() < 1e-12: continue
            u3 = G3.move(zm, Cm, G3.W_(h, y13))
            for gi, g in enumerate(MR):
                if abs(cs[gi]) > 1e-12: lhs += cs[gi] * fl(G3, G3.W_(g, u3))
        c, *_ = np.linalg.lstsq(B2, lhs, rcond=None)
        rel = np.linalg.norm(B2 @ c - lhs) / max(np.linalg.norm(lhs), 1e-300)
        big = np.argsort(-np.abs(x))[:8]
        terms = [(MR[i % len(MR)], MC[i // len(MR)], round(float(abs(x[i])), 3)) for i in big if abs(x[i]) > 1e-6]
        print(f"  identity {k}: at f={f3} distance from L2 {rel:.2e}  terms (g, h, |c|): {terms}", flush=True)


if __name__ == "__main__":
    main()
