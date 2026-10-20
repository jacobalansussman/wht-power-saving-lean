"""s6_core_probe.py -- block-4 cores in the no-scratch reflection form (evidence only, calibrated).

A core is a no-scratch reflection program at block 4 (15 gates, each direction once) whose targets R_l are any
symmetric matrices with sum 0 (det condition) and total plain length sum len(R_l) >= 19 > 15, i.e. a block-4
saving program with mixed/alternating targets. (With phase letters at the ends, alternating nonsingular targets are
engine targets, rm2-phase Lemma A; four of them in 15 would be the Reed-Muller member 4 in 15, saving 0.047.)

Modes:  altpos  -- 4 distinct alternating targets, plain 20 gates (each target on one of its 5-direction frames):
                   positive control, must be found;
        alt4    -- 4 distinct alternating targets summing to 0, 15 gates (every direction once), random orders;
        mix     -- targets 3 invertible non-alternating + 1 alternating + 1 rank-2 (len 4+4+4+5+2 = 19), sum 0.
"""
import sys, itertools, time
import numpy as np
sys.path.insert(0, ".")
from s4_refl_probe import setup, chars, lm

m = 4
dirs, zbits, pairs, Qs, QZ = setup(m)

def f2rank(T):
    Mm = T.copy() % 2; r = 0
    for c in range(Mm.shape[1]):
        p = next((i for i in range(r, Mm.shape[0]) if Mm[i, c]), None)
        if p is None: continue
        Mm[[r, p]] = Mm[[p, r]]
        for i in range(Mm.shape[0]):
            if i != r and Mm[i, c]: Mm[i] ^= Mm[r]
        r += 1
    return r

allsym = []
for bits in itertools.product((0, 1), repeat=10):
    T = np.zeros((4, 4), dtype=np.int64)
    for b, (i, j) in zip(bits, pairs): T[i, j] = T[j, i] = b
    allsym.append(T)
alt = [T for T in allsym if f2rank(T) == 4 and not T.diagonal().any()]
inv_na = [T for T in allsym if f2rank(T) == 4 and T.diagonal().any()]
rank2 = [T for T in allsym if f2rank(T) == 2]
print("alternating nonsingular:", len(alt), " invertible non-alternating:", len(inv_na), " rank 2:", len(rank2))

def frames_of(T, maxn=5):
    # all sets of rank-one matrices (directions) of size len(T) summing to T (search up to size maxn)
    out = []
    for k in range(1, maxn + 1):
        for comb in itertools.combinations(dirs, k):
            S = sum(np.outer(zbits[z - 1], zbits[z - 1]) for z in comb) % 2
            if (S == T % 2).all(): out.append(list(comb))
        if out: return out
    return out

def run(mode, seed, nstarts=4):
    rng = np.random.default_rng(seed)
    if mode == "altpos":
        idx = rng.choice(len(alt), size=4, replace=False); Ts = [alt[i] for i in idx]
        gl = []
        for T in Ts:
            fr = frames_of(T); gl += fr[rng.integers(len(fr))]
        gl = [gl[i] for i in rng.permutation(len(gl))]
    elif mode == "alt4":
        while True:
            idx = rng.choice(len(alt), size=3, replace=False); Ts = [alt[i] for i in idx]
            T4 = sum(Ts) % 2
            if f2rank(T4) == 4 and not T4.diagonal().any() and all((T4 != t).any() for t in Ts):
                Ts.append(T4); break
        gl = list(rng.permutation(dirs))
    elif mode == "mix":
        while True:
            Ts = [inv_na[i] for i in rng.choice(len(inv_na), size=3, replace=False)] + [alt[rng.integers(len(alt))]]
            T5 = sum(Ts) % 2
            if f2rank(T5) == 2: Ts.append(T5); break
        gl = list(rng.permutation(dirs))
    W = len(Ts)
    lens = [len(frames_of(T)[0]) for T in Ts]
    act = np.stack([QZ[:, z - 1] for z in gl], axis=1)
    D = np.zeros((len(Qs), W, W))
    for l, T in enumerate(Ts): D[:, l, l] = 1 - 2 * chars(Qs, pairs, T)
    print("%s seed %d: W = %d registers, %d gates, target lengths %s (sum %d)" % (mode, seed, W, len(gl), lens, sum(lens)), flush=True)
    for start in range(nstarts):
        U, f, it = lm(rng.normal(size=(len(gl), W)), act, D)
        print("   start %d: residual %.6g after %d its" % (start, f, it), flush=True)

if __name__ == "__main__":
    mode = sys.argv[1]; seed = int(sys.argv[2]); n = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    t0 = time.time(); run(mode, seed, n); print("time %.1f s" % (time.time() - t0))
