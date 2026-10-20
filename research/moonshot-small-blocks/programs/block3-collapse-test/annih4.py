"""annih4.py -- ms-xw3.  Layer-collapse test with mixes at BOTH internal positions of a 3-move word:
columns  g . M3 . h . M2 . k . M1 (v)   (g in R, h in H, k in Kset; M1, M2, M3 fixed paid moves, possibly in frames)
against  L2 = span{ b(x) : b in B, x in {v, M_i v, M_j k M_i v (i != j, k in Kset)} }  (every word with <= 2 of the
moves and a mix from Kset between them).  f1 and f2 stacked, one coefficient vector; EXCESS = degenerate - trivial.
usage: python3 -I annih4.py <moves z,C;z,C;z,C[|z,C|z,C ...]> <R list> <H list|menu> <K list> <B list|menu> <f1,f2> <S>
Extra top moves after '|' are added as further third moves (several top moves at once: the shared-middle shape, where
the other array's work reaches an output through two different private top moves)."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw3lib import Game, menu as MENU
from annih2 import fl


def L(arg):
    try: return MENU(arg)
    except ValueError: return arg.split(",")


def deg_nontrivial(A, Ap, triv):
    """dimension of span(non-trivial columns) cap L2, from the full and projected column sets"""
    import numpy as _np
    keep = _np.setdiff1d(_np.arange(A.shape[1]), triv)
    if keep.size == 0: return 0
    s1 = _np.linalg.svd(A[:, keep], compute_uv=False); r1 = int((s1 > s1[0] * 1e-9).sum())
    s2 = _np.linalg.svd(Ap[:, keep], compute_uv=False); r2 = int((s2 > s2[0] * 1e-9).sum())
    return r1 - r2


def main():
    parts = sys.argv[1].split("|")
    mv = [((m.split(",") + ["1"])[0], (m.split(",") + ["1"])[1]) for m in parts[0].split(";")]
    tops = [mv[2]] + [((m.split(",") + ["1"])[0], (m.split(",") + ["1"])[1]) for m in parts[1:]]
    mv = mv + tops[1:]
    R, H, Ks, B = L(sys.argv[2]), L(sys.argv[3]), L(sys.argv[4]), L(sys.argv[5])
    FS = [int(a) for a in sys.argv[6].split(",")]; S = int(sys.argv[7])
    t0 = time.time()
    G = Game(FS, W=1, S=S, seed=11)
    v = G.inp(0)
    M = lambda i, X: G.move(int(mv[i][0]), mv[i][1], X)
    y1 = M(0, v)
    cols = []; lab = []
    tidx = [2] + list(range(3, len(mv)))
    for t in tidx:
        for k in Ks:
            y2 = M(1, G.W_(k, y1))
            for h in H:
                y3 = M(t, G.W_(h, y2))
                for g in R:
                    cols.append(fl(G, G.W_(g, y3))); lab.append((t, g, h, k))
    A3 = np.stack(cols, 1)
    nm = len(mv)
    xs = [v] + [M(i, v) for i in range(nm)]
    for i in range(nm):
        for j in range(nm):
            if i != j:
                for k in Ks: xs.append(M(j, G.W_(k, M(i, v))))
    A2 = np.stack([fl(G, G.W_(b, x)) for x in xs for b in B], 1)
    U2, s2, _ = np.linalg.svd(A2, full_matrices=False); Q2 = U2[:, s2 > s2[0] * 1e-9]
    del U2
    A3p = A3 - Q2 @ (Q2.conj().T @ A3)
    sA = np.linalg.svd(A3, compute_uv=False); rA = int((sA > sA[0] * 1e-9).sum())
    sP = np.linalg.svd(A3p, compute_uv=False); rP = int((sP > sP[0] * 1e-9).sum())
    cn = np.linalg.norm(A3p, axis=0) / np.maximum(np.linalg.norm(A3, axis=0), 1e-300)
    triv = np.flatnonzero(cn < 1e-8); rT = 0
    if triv.size:
        sT = np.linalg.svd(A3[:, triv], compute_uv=False); rT = int((sT > sT[0] * 1e-9).sum())
    print(f"# annih4 moves={sys.argv[1]} |R|={len(R)} |H|={len(H)} |K|={len(Ks)} |B|={len(B)} FS={FS} S={S}: cols {A3.shape[1]}, "
          f"L2 cols {A2.shape[1]} (rank {Q2.shape[1]}), rows {A3.shape[0]}; rank A3 {rA}, mod L2 {rP}: degenerate {rA - rP}, "
          f"trivial {rT}, EXCESS {rA - rP - rT}, NONTRIVIAL-COLLAPSE {deg_nontrivial(A3, A3p, triv)}; gap {sP[rP-1]/sP[0]:.1e} / {sP[rP]/sP[0] if rP < sP.size else 0:.1e} "
          f"({time.time()-t0:.0f}s)", flush=True)
    if triv.size:
        print("#   trivial (g,h,k) examples:", [lab[i] for i in triv[:8]], flush=True)


if __name__ == "__main__":
    main()
