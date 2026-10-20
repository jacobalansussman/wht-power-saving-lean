"""v-ms-xw3: collapse test for an n-move path with mixes, against the COMPLETE span of every term with fewer moves.

Path: v -> M1 -> (mix W1) -> M2 -> (mix W2) -> ... -> Mn -> (readout R).  The content of move j+1 is a mix from W_j
of ANY earlier result, so the lower terms are all sub-paths (subsequences of the moves, at least one move dropped):
 - if Mn is kept, the term is  R Mn W.. M.. (readout R in front);
 - if not, the output reads the last kept move directly through B;
 - the first kept move M_i reads v through W_{i-1} (M1 reads v through the identity, as in ms-xw3's tests).
This is the full lower span of ms-xw3's annih4 (n = 3) and annih6 (n = 4) shapes; theirs omits the sub-paths through
the top move that skip a lower move (e.g. R M3 H M1 v, R M3 H v, R M3 H M2 K v).
usage: path_ext.py <moves z,C;z,C;...> <R> <W1> ... <W_{n-1}> <B> <f1,f2> <S> [f3]
"""
import sys, os, time, itertools
import numpy as np
np.seterr(all="ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xtalk_ext import L, fl, zc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xw3lib import Game


def path_vecs(G, v, mv, W, chosen, front):
    """all vectors front . M_{c_k} . W_{c_k - 1} . ... . M_{c_1} . W_{c_1 - 1} . v  (W_0 = identity), with labels"""
    cur = [(v, ())]
    for idx, i in enumerate(chosen):
        mixes = ["1"] if i == 0 else W[i - 1]
        nxt = []
        for (x, lab) in cur:
            for k in mixes:
                nxt.append((G.move(mv[i][0], mv[i][1], G.W_(k, x)), lab + (k,)))
        cur = nxt
    out = []
    for (x, lab) in cur:
        for g in front:
            out.append((fl(G, G.W_(g, x)), lab + (g,)))
    return out


def build(G, mv, R, W, B):
    n = len(mv); v = G.inp(0)
    top = path_vecs(G, v, mv, W, list(range(n)), R)
    low = []
    for k in range(n):
        for sub in itertools.combinations(range(n), k):
            front = R if (k > 0 and sub[-1] == n - 1) else B
            low += path_vecs(G, v, mv, W, list(sub), front)
    return top, low


def orth(A):
    U, s, _ = np.linalg.svd(A, full_matrices=False); r = int((s > s[0] * 1e-9).sum()); return U[:, :r]


def rank(M):
    if M.shape[1] == 0: return 0, np.array([1.0])
    s = np.linalg.svd(M, compute_uv=False); return int((s > s[0] * 1e-9).sum()), s


def main():
    mv = [zc(m) for m in sys.argv[1].split(";")]; n = len(mv)
    R = L(sys.argv[2]); W = [L(sys.argv[3 + j]) for j in range(n - 1)]; B = L(sys.argv[2 + n])
    FS = [int(x) for x in sys.argv[3 + n].split(",")]; S = int(sys.argv[4 + n])
    f3 = int(sys.argv[5 + n]) if len(sys.argv) > 5 + n else 5
    t0 = time.time()
    G = Game(FS, W=1, S=S, seed=31)
    top, low = build(G, mv, R, W, B)
    A3 = np.stack([t[0] for t in top], 1); lab = [t[1] for t in top]
    A2 = np.stack([t[0] for t in low], 1); del low
    print(f"# path_ext moves={sys.argv[1]} |R|={len(R)} |W|={[len(w) for w in W]} |B|={len(B)} FS={FS} S={S}: rows {A3.shape[0]}, "
          f"top cols {A3.shape[1]}, lower cols {A2.shape[1]}", flush=True)
    if A3.shape[1] + A2.shape[1] > 0.9 * A3.shape[0]: print("  WARNING: columns close to rows", flush=True)
    Q = orth(A2); del A2
    A3p = A3 - Q @ (Q.conj().T @ A3)
    rA, sA = rank(A3); rP, sP = rank(A3p)
    cn = np.linalg.norm(A3p, axis=0) / np.linalg.norm(A3, axis=0)
    triv = np.flatnonzero(cn < 1e-8); keep = np.setdiff1d(np.arange(A3.shape[1]), triv)
    rT = rank(A3[:, triv])[0] if triv.size else 0
    nt = (rank(A3[:, keep])[0] - rank(A3p[:, keep])[0]) if keep.size else 0
    if nt:
        _, s2, Vh = np.linalg.svd(A3p[:, keep], full_matrices=True)
        r2 = int((s2 > s2[0] * 1e-9).sum()); N = Vh[r2:].conj().T
        Wm = A3[:, keep] @ N
        u, sw, vh = np.linalg.svd(Wm, full_matrices=False)
        nt_ex = int((sw > sA[0] * 1e-8).sum()); print(f"  extracted {nt_ex}", flush=True); nt = min(nt, nt_ex)
    print(f"  lower rank {Q.shape[1]}; rank top {rA}, mod lower {rP}: degenerate {rA - rP}, trivial {rT}, NONTRIVIAL {nt}; "
          f"sv gap {sP[rP-1]/sP[0]:.1e} / {(sP[rP]/sP[0]) if rP < sP.size else 0:.1e} ({time.time()-t0:.0f}s)", flush=True)
    if triv.size: print("  trivial examples:", [lab[i] for i in triv[:6]], flush=True)
    if nt:
        Xg = N @ vh[:nt].conj().T
        G3 = Game([f3], W=1, S=1, seed=78)
        top3, low3 = build(G3, mv, R, W, B)
        B3 = np.stack([t[0] for t in top3], 1)[:, keep]; Q3 = orth(np.stack([t[0] for t in low3], 1))
        for j in range(min(nt, 6)):
            x = Xg[:, j]; lhs = B3 @ x; res = lhs - Q3 @ (Q3.conj().T @ lhs)
            big = np.argsort(-np.abs(x))[:10]
            print(f"   collapse {j}: at f={f3} distance from complete lower span {np.linalg.norm(res)/np.linalg.norm(lhs):.2e}; terms "
                  f"{[(lab[keep[i]], np.round(abs(x[i]), 3)) for i in big if abs(x[i]) > 1e-6]}", flush=True)


if __name__ == "__main__":
    main()
