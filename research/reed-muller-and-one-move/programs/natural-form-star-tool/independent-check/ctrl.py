"""ctrl.py: genuine plain programs (canonical frames, disjoint) at block 5 with one move per direction; build the TRUE
E = I + Gamma - K of the program (Lemma 1 with identity mixings), check T10 on all 32,767 forms, then run vrec.run
with trueE to check that every bound u(A) >= true rank and every derived zero cell is a true zero."""
import sys, itertools, numpy as np
import vf, vrec
def frame(k): return [vf.coord(vf.gmul(vf.gpow(vf.ALPHA, k), f)) for f in vf.BASIS]
def trueE(ks, order_seed=0):
    rng = np.random.default_rng(order_seed)
    E = np.eye(31)
    for k in ks:
        fr = list(frame(k)); rng.shuffle(fr)          # time order inside the frame
        for a, x in enumerate(fr):
            E[x - 1, x - 1] = 0
            for b, y in enumerate(fr):
                if b < a: E[x - 1, y - 1] = 1          # x later than y: lower entry
                if b > a: E[x - 1, y - 1] = -1         # x earlier: upper entry
    return E
def check_T10(ks, E):
    w = vf.PAR[:, ks].sum(axis=1)
    bad = 0
    for Q in range(1, vf.NF):
        D = [z - 1 for z in range(1, 32) if int(vf.MASK[Q]) >> (z - 1) & 1]
        rk = np.linalg.matrix_rank(E[np.ix_(D, D)])
        if rk != len(D) - w[Q]: bad += 1
    return bad
if __name__ == "__main__":
    sets = [[0], [0, 3], [0, 3, 6], [2, 3, 20, 30], [4, 11, 15, 24, 29]]
    for ks in sets:
        fr = [frame(k) for k in ks]; allz = [z for f in fr for z in f]
        assert len(set(allz)) == len(allz), ("frames overlap", ks)
        E = trueE(ks)
        print(ks, "T10 mismatches:", check_T10(ks, E), flush=True)
        print(ks, "vrec with true E:", vrec.run(ks, maxd=int(sys.argv[1]) if len(sys.argv) > 1 else 12, rounds=3, trueE=E), flush=True)
