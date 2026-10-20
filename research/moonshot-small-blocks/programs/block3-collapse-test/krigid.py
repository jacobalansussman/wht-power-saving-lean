"""krigid.py -- ms-xw3.  Zero-divisor test for readouts: is there a nonzero mix rho of free words with
rho . K  in  L2 = span{ b(y) : b in menuB, y a <= 2-layer product of D1, D2, D4 (and of frame moves) applied to v }?
If not, no readout of a shared top move can be a zero divisor (left annihilator argument, NOTE section 3).
usage: python3 -I krigid.py <menuR> <menuB> <f1,f2> <S> [extra frame moves 'z,C;z,C']"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from xw3lib import Game, menu as MENU
from annih2 import fl


def main():
    MR, MB = MENU(sys.argv[1]), MENU(sys.argv[2]); FS = [int(a) for a in sys.argv[3].split(",")]; S = int(sys.argv[4])
    extra = [tuple((m.split(",") + ["1"])[:2]) for m in sys.argv[5].split(";")] if len(sys.argv) > 5 else []
    t0 = time.time()
    G = Game(FS, W=1, S=S, seed=9)
    v = G.inp(0)
    Kv = G.K(v)
    mv = [("1", "1"), ("2", "1"), ("4", "1")] + extra
    ys = [v] + [G.move(int(z), C, v) for z, C in mv]
    ys += [G.move(int(z2), C2, G.move(int(z1), C1, v)) for i, (z1, C1) in enumerate(mv) for (z2, C2) in mv[i + 1:]]
    A = np.stack([fl(G, G.W_(g, Kv)) for g in MR], 1)
    B = np.stack([fl(G, G.W_(b, y)) for y in ys for b in MB], 1)
    U, s, _ = np.linalg.svd(B, full_matrices=False); Q = U[:, s > s[0] * 1e-9]
    Ap = A - Q @ (Q.conj().T @ A)
    sA = np.linalg.svd(A, compute_uv=False); rA = int((sA > sA[0] * 1e-9).sum())
    sP = np.linalg.svd(Ap, compute_uv=False); rP = int((sP > sP[0] * 1e-9).sum())
    print(f"# krigid R={sys.argv[1]}({len(MR)}) B={sys.argv[2]}({len(MB)}) <=2-layer sources {len(ys)} FS={FS} S={S} rows {A.shape[0]} "
          f"L2 cols {B.shape[1]} (rank {Q.shape[1]}): rank(rho K) {rA}, mod L2 {rP} -> zero-divisor dimension {rA - rP}; "
          f"smallest kept sv {sP[rP-1]/sP[0]:.2e}, next {sP[rP]/sP[0] if rP < sP.size else 0:.2e} ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
