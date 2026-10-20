"""v-ms-xw3: the shared-top cross-talk test of ms-xw3 (annih2/annih3), redone with the COMPLETE lower span.

ms-xw3's L2 = span{ b(y1), b(y0), b(v) : b in B }.  In the shared-top shape the shared move m also reads the other
array's LOWER results, so the output also contains  g m h(y0)  (two moves) and  g m h(v)  (one move), with g a readout
word and h a content word.  Those are terms with fewer paid moves too, and a collapse may use them.  Here
L2ext = L2 + span{ g m_k h(y0), g m_k h(v) : g in R, h in C, every top move m_k }.
Reported for BOTH L2 and L2ext: degenerate = dim(span A3 cap L), trivial = rank of columns individually in L,
nontrivial = dim(span(other columns) cap L).
usage: xtalk_ext.py <R list|menu> <C list|menu> <B list|menu> <chain> <moves 'z,C;z,C'> <f1,f2> <S>"""
import sys, os, time
import numpy as np
np.seterr(all="ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))  # fallback copy (Colab)
from xw3lib import Game, menu as MENU


def L(arg):
    try: return MENU(arg)
    except ValueError: return arg.split(",")


def fl(G, X):
    return np.concatenate([X[f].reshape(-1) for f in G.FS])


def zc(s):
    p = (s.split(",") + ["1"])[:2]; return int(p[0]), p[1]


def chain_vecs(G, chain):
    v = G.inp(0)
    q = chain.split(";")
    if len(q) == 1:
        z0, C0 = zc(q[0]); return v, v, G.move(z0, C0, v)
    z0, C0 = zc(q[0]); z1, C1 = zc(q[2])
    y0 = G.move(z0, C0, v)
    return v, y0, G.move(z1, C1, G.W_(q[1], y0))


def rank(M, tol=1e-9):
    if M.shape[1] == 0: return 0, np.array([])
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > s[0] * tol).sum()), s


def analyse(A3, A2, tag):
    U, s, _ = np.linalg.svd(A2, full_matrices=False)
    r2 = int((s > s[0] * 1e-9).sum()); Q = U[:, :r2]; del U
    A3p = A3 - Q @ (Q.conj().T @ A3)
    rA, _ = rank(A3); rP, sP = rank(A3p)
    cn = np.linalg.norm(A3p, axis=0) / np.maximum(np.linalg.norm(A3, axis=0), 1e-300)
    triv = np.flatnonzero(cn < 1e-8)
    rT = rank(A3[:, triv])[0] if triv.size else 0
    keep = np.setdiff1d(np.arange(A3.shape[1]), triv)
    nt = (rank(A3[:, keep])[0] - rank(A3p[:, keep])[0]) if keep.size else 0
    gap = (sP[rP - 1] / sP[0], sP[rP] / sP[0] if rP < sP.size else 0.0)
    print(f"  [{tag}] L cols {A2.shape[1]} rank {r2}; rank A3 {rA}, mod L {rP}: degenerate {rA - rP}, trivial {rT}, "
          f"NONTRIVIAL {nt}; sv gap {gap[0]:.1e} / {gap[1]:.1e}", flush=True)
    return rA - rP, rT, nt


def main():
    R, C, B = L(sys.argv[1]), L(sys.argv[2]), L(sys.argv[3]); chain = sys.argv[4]
    moves = [zc(m) for m in sys.argv[5].split(";")]
    FS = [int(x) for x in sys.argv[6].split(",")]; S = int(sys.argv[7])
    t0 = time.time()
    G = Game(FS, W=1, S=S, seed=21)
    v, y0, y1 = chain_vecs(G, chain)
    cols = []
    for (z, Cm) in moves:
        for h in C:
            u = G.move(z, Cm, G.W_(h, y1))
            for g in R: cols.append(fl(G, G.W_(g, u)))
    A3 = np.stack(cols, 1)
    base = [fl(G, G.W_(b, s)) for s in ((y1, y0, v) if y0 is not v else (y1, v)) for b in B]
    ext = []
    for (z, Cm) in moves:
        for h in C:
            for s in ((y0, v) if y0 is not v else (v,)):
                u = G.move(z, Cm, G.W_(h, s))
                for g in R: ext.append(fl(G, G.W_(g, u)))
    A2 = np.stack(base, 1); A2x = np.stack(base + ext, 1)
    print(f"# xtalk_ext R={sys.argv[1]}({len(R)}) C={sys.argv[2]}({len(C)}) B={sys.argv[3]}({len(B)}) chain={chain} "
          f"moves={sys.argv[5]} FS={FS} S={S} rows {A3.shape[0]} A3 cols {A3.shape[1]} ext cols {len(ext)}", flush=True)
    if A3.shape[1] + A2x.shape[1] > 0.9 * A3.shape[0]:
        print("  WARNING: columns close to rows; raise S", flush=True)
    analyse(A3, A2, "ms-xw3 L2")
    analyse(A3, A2x, "complete L2ext")
    print(f"  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
