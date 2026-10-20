"""annih6.py -- ms-xw3.  One level up (two shared moves in sequence, the 3-in-8 shape where the second shared move
reads the first): columns  g . t1 . h . t0 . k (Y),  Y = M2 M1 v a private 2-move chain, against
L3 = span{ b(x) : x in {v, M1 v, Y, t(k Y), t(k M1 v), t(k v) : t in {t0, t1}, k in Kset} }  (every <= 3-layer term such a
unit has).  EXCESS = degenerate - trivial, as in annih2/4.
usage: python3 -I annih6.py <chain z,C;z,C> <t0 z,C> <t1 z,C> <R list> <H list> <K list> <B list> <f1,f2> <S>"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw3lib import Game, menu as MENU
from annih2 import fl


def L(arg):
    try: return MENU(arg)
    except ValueError: return arg.split(",")


def zc(s):
    a = (s.split(",") + ["1"])[:2]; return int(a[0]), a[1]


def deg_nontrivial(A, Ap, triv):
    """dimension of span(non-trivial columns) cap L2, from the full and projected column sets"""
    import numpy as _np
    keep = _np.setdiff1d(_np.arange(A.shape[1]), triv)
    if keep.size == 0: return 0
    s1 = _np.linalg.svd(A[:, keep], compute_uv=False); r1 = int((s1 > s1[0] * 1e-9).sum())
    s2 = _np.linalg.svd(Ap[:, keep], compute_uv=False); r2 = int((s2 > s2[0] * 1e-9).sum())
    return r1 - r2


def main():
    ch = [zc(q) for q in sys.argv[1].split(";")]; t0, t1 = zc(sys.argv[2]), zc(sys.argv[3])
    R, H, Ks, B = L(sys.argv[4]), L(sys.argv[5]), L(sys.argv[6]), L(sys.argv[7])
    FS = [int(a) for a in sys.argv[8].split(",")]; S = int(sys.argv[9])
    tm = time.time()
    G = Game(FS, W=1, S=S, seed=13)
    v = G.inp(0)
    y1 = G.move(ch[0][0], ch[0][1], v); Y = G.move(ch[1][0], ch[1][1], y1)
    cols = []; lab = []
    for k in Ks:
        u0 = G.move(t0[0], t0[1], G.W_(k, Y))
        for h in H:
            u1 = G.move(t1[0], t1[1], G.W_(h, u0))
            for g in R:
                cols.append(fl(G, G.W_(g, u1))); lab.append((g, h, k))
    A4 = np.stack(cols, 1)
    xs = [v, y1, Y]
    for t in (t0, t1):
        for k in Ks:
            for src in (Y, y1, v):
                xs.append(G.move(t[0], t[1], G.W_(k, src)))
    A3 = np.stack([fl(G, G.W_(b, x)) for x in xs for b in B], 1)
    U2, s2, _ = np.linalg.svd(A3, full_matrices=False); Q = U2[:, s2 > s2[0] * 1e-9]; del U2
    A4p = A4 - Q @ (Q.conj().T @ A4)
    sA = np.linalg.svd(A4, compute_uv=False); rA = int((sA > sA[0] * 1e-9).sum())
    sP = np.linalg.svd(A4p, compute_uv=False); rP = int((sP > sP[0] * 1e-9).sum())
    cn = np.linalg.norm(A4p, axis=0) / np.maximum(np.linalg.norm(A4, axis=0), 1e-300)
    triv = np.flatnonzero(cn < 1e-8); rT = 0
    if triv.size:
        sT = np.linalg.svd(A4[:, triv], compute_uv=False); rT = int((sT > sT[0] * 1e-9).sum())
    print(f"# annih6 chain={sys.argv[1]} t0={sys.argv[2]} t1={sys.argv[3]} |R|={len(R)} |H|={len(H)} |K|={len(Ks)} |B|={len(B)} "
          f"FS={FS} S={S}: cols {A4.shape[1]}, L3 cols {A3.shape[1]} (rank {Q.shape[1]}), rows {A4.shape[0]}; rank {rA}, "
          f"mod L3 {rP}: degenerate {rA - rP}, trivial {rT}, EXCESS {rA - rP - rT}, NONTRIVIAL-COLLAPSE {deg_nontrivial(A4, A4p, triv)}; gap {sP[rP-1]/sP[0]:.1e} / "
          f"{sP[rP]/sP[0] if rP < sP.size else 0:.1e} ({time.time()-tm:.0f}s)", flush=True)
    if triv.size: print("#   trivial (g,h,k) examples:", [lab[i] for i in triv[:8]], flush=True)


if __name__ == "__main__":
    main()
