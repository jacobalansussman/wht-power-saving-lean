"""vscan.py (v-r3-3in14): my own scan of plain 3 in 15 minus one move at block 5, in a basis-free model
(field elements as directions, the polynomial basis for Sym).  Symmetry: e -> 2^f e + j (mod 31) on the logs of
directions AND on target labels k (c = alpha^(2k)); derivation in NOTE.md section 1.
Output: out_vscan.txt, vsets.json (the passing sets with my zero-block verdicts)."""
import itertools, json, sys
import numpy as np
from vcore import EXP, LOG, dir_vec, tgt_vec, parity_table, Inst

DIRS = list(range(1, 32))                       # nonzero field elements (polynomial-basis ints)
DV = {x: dir_vec(x) for x in DIRS}
TV = {k: tgt_vec(EXP[(2 * k) % 31]) for k in range(31)}
DT = {x: parity_table(DV[x]) for x in DIRS}
TT = {k: parity_table(TV[k]) for k in range(31)}

# sanity: the 31 direction forms span all of Sym (dimension 15)
def rank2(vs):
    vs = list(vs); r = 0
    for b in range(15):
        piv = next((i for i in range(r, len(vs)) if vs[i] >> b & 1), None)
        if piv is None: continue
        vs[r], vs[piv] = vs[piv], vs[r]
        for i in range(len(vs)):
            if i != r and vs[i] >> b & 1: vs[i] ^= vs[r]
        r += 1
    return r
assert rank2(DV.values()) == 15

# frames: 5-sets of distinct directions whose forms add up to the target's form
frames = {k: [] for k in range(31)}
inv = {v: k for k, v in TV.items()}
assert len(inv) == 31
for S in itertools.combinations(DIRS, 5):
    v = 0
    for x in S: v ^= DV[x]
    if v in inv: frames[inv[v]].append(S)
nfr = sorted(set(len(f) for f in frames.values()))

G = [(f, j) for f in range(5) for j in range(31)]
def act(g, e):
    f, j = g
    return ((e << f) + j) % 31
def canon(ks, logs):
    best = None
    for g in G:
        key = (tuple(sorted(act(g, k) for k in ks)), tuple(sorted(act(g, e) for e in logs)))
        if best is None or key < best: best = key
    return best

ms = sorted(set(canon(c, []) [0] for c in itertools.combinations_with_replacement(range(31), 3)))
reps = ms
plain_keys = set(); passing = {}
for ks in reps:
    for fr in itertools.product(*(range(len(frames[k])) for k in ks)):
        moves = [x for k, f in zip(ks, fr) for x in frames[k][f]]
        cnt = sum(DT[x].astype(np.int32) for x in moves)
        r = sum(TT[k].astype(np.int32) for k in ks)
        slack = cnt - r
        assert (slack >= 0).all(), "plain fails deletion?!"
        logs = [LOG[x] for x in moves]
        plain_keys.add(canon(ks, logs))
        tight = slack == 0
        for d in sorted(set(moves)):
            if (tight & (DT[d] == 1)).any(): continue          # deleting d breaks deletion
            m14 = list(moves); m14.remove(d)
            key = canon(ks, [LOG[x] for x in m14])
            if key not in passing:
                passing[key] = {"ks": list(ks), "frames": list(fr), "delete": d, "mult": moves.count(d),
                                "moves": m14, "tight_before": int(tight.sum())}

out = []
out.append(f"frames per target: {nfr}")
out.append(f"classes of 3-multisets of targets: {len(reps)}")
out.append(f"plain multisets (classes): {len(plain_keys)}")
out.append(f"plain minus one passing span+deletion (classes): {len(passing)}")
res = []
nex = 0
for key, s in sorted(passing.items()):
    I = Inst([TV[k] for k in s["ks"]], [DV[x] for x in s["moves"]])
    assert I.deletion_ok()
    ex, cert, c = I.zero_block()
    nex += ex
    s.update({"key": [list(key[0]), list(key[1])], "excluded": ex, "cert": cert, "clique": c})
    res.append(s)
    out.append(f"ks={tuple(s['ks'])} mult={s['mult']} tight_before={s['tight_before']} key={key} clique={c} "
               f"Z={'EXCLUDED' if ex else 'pass'} cert={cert and {k: cert[k] for k in ('kind','F','DF','rF','Y') if k in cert}}")
out.append(f"zero-block: {nex} of {len(res)} excluded")
open("out_vscan.txt", "w").write("\n".join(out) + "\n")
json.dump(res, open("vsets.json", "w"))
print("\n".join(out[:4] + out[-1:]))
