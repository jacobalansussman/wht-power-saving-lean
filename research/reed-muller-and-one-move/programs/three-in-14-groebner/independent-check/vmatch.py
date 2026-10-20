"""vmatch.py (v-r3-3in14): (1) convert the track's instances (its certs_min.json: its own 44 in basis
[3,5,12,17,26], rm2-solve's 44 in rm2-solve's coordinates, the controls) to field elements, by finding for each
coordinate system the ordered self-dual basis that reproduces ITS target matrices exactly; (2) compare the 44 + 44
with my own 44 class by class (canonical key under e -> 2^f e + j); (3) run MY zero-block test on every one of
their instances in THEIR coordinates (generic code on their matrices), and re-check THEIR certificates with my own
checker.  Output: out_vmatch.txt."""
import itertools, json
from vcore import EXP, LOG, mul, tr, Inst, mat_to_vec, zz_vec

TRK = "../"
D = json.load(open(TRK + "certs_min.json"))
mine = json.load(open("vsets.json"))
mykeys = {tuple(map(tuple, s["key"])) for s in mine}

# all ordered self-dual bases of GF(32)
elems = list(range(1, 32))
ON = [x for x in elems if tr(mul(x, x)) == 1]
SDB = []
def rec(cur):
    if len(cur) == 5:
        SDB.append(tuple(cur)); return
    for x in ON:
        if x in cur: continue
        if all(tr(mul(x, y)) == 0 for y in cur): rec(cur + [x])
rec([])

def Tmat(c, f):
    return [[tr(mul(c, mul(f[i], f[j]))) for j in range(5)] for i in range(5)]

def find_basis(ks, Ts):
    hits = [f for f in SDB if all(Tmat(EXP[(2 * k) % 31], f) == T for k, T in zip(ks, Ts))]
    return hits

G = [(fq, j) for fq in range(5) for j in range(31)]
def canon(ks, logs):
    return min((tuple(sorted(((k << fq) + j) % 31 for k in ks)), tuple(sorted(((e << fq) + j) % 31 for e in logs)))
               for fq, j in G)

def to_field(z, f):
    x = 0
    for j in range(5):
        if z >> j & 1: x ^= f[j]
    return x

out = [f"ordered self-dual bases: {len(SDB)}"]
bykey = {}
nZ = {}; ncert = {}
for inst in D:
    tag = inst["tag"]; m = inst["m"]
    if m != 5:
        continue
    hits = find_basis(inst["ks"], inst["Ts"])
    # the basis must also be consistent with the frames/moves (any hit reproduces the targets; the targets with
    # k = 0 and others fix the basis only up to the stabiliser, so try all hits and keep the canonical keys)
    keys = set()
    for f in hits:
        xs = [to_field(z, f) for z in inst["moves"]]
        assert all(xs)
        keys.add(canon(inst["ks"], [LOG[x] for x in xs]))
    bykey.setdefault(tag, []).append((inst["idx"], keys))
    # my zero-block test on THEIR coordinates and matrices
    I = Inst([mat_to_vec(T) for T in inst["Ts"]], [zz_vec(z) for z in inst["moves"]])
    ok = I.deletion_ok()
    ex, cert, c = I.zero_block()
    verdict = "EXCLUDED" if (ex or not ok) else "passes"
    nZ.setdefault(tag, [0, 0]); nZ[tag][0] += (verdict == inst["verdict"]); nZ[tag][1] += 1
    # check THEIR certificate with my own checker (Z kind)
    tc = inst.get("cert")
    if tc and tc.get("kind") == "Z":
        Y = tc["Y"]; F = tc["F"]
        # their F is a form in their own monomial encoding; recompute D_F from MY side via their DF mask, and check
        # (a) every pair in Y lies in a tight support of MY tight list; (b) |D_F| + r_F < 2|Y| with D_F, r_F = MY
        # values for some form whose support is exactly their DF.
        ts = I.tight_supports()
        pairs_ok = all(any(t >> i & 1 and t >> j & 1 for t in ts) for i in Y for j in Y)
        DF = tc["DF"]
        import numpy as np
        qs = np.nonzero(I.Dmask == DF)[0]
        rs = sorted(set(int(I.r[q]) for q in qs))
        dsz = bin(DF).count("1")
        ok2 = pairs_ok and len(qs) > 0 and all(Y_i >= 0 for Y_i in Y) and (dsz + min(rs) < 2 * len(Y)) and \
              all(DF >> i & 1 for i in Y)
        ncert.setdefault(tag, [0, 0]); ncert[tag][0] += ok2; ncert[tag][1] += 1
        if inst["idx"] < 3 or not ok2:
            out.append(f"  {tag} idx {inst['idx']}: their cert Y={Y} |DF|={dsz} r(forms with that support)={rs} "
                       f"pairs-in-my-tight-supports={pairs_ok} -> {'OK' if ok2 else 'BAD'}")
for tag, v in nZ.items():
    out.append(f"{tag}: my zero-block verdict agrees with theirs on {v[0]} of {v[1]}")
for tag, v in ncert.items():
    out.append(f"{tag}: their Z certificate re-checked by my code: {v[0]} of {v[1]} valid")
for tag in ("mine-b5", "rm2-b5"):
    ks_ = [k for _, kk in bykey[tag] for k in kk]
    amb = sum(1 for _, kk in bykey[tag] if len(kk) != 1)
    S = set(ks_)
    out.append(f"{tag}: {len(bykey[tag])} instances, {amb} with an ambiguous basis, {len(S)} distinct classes; "
               f"equal to my 44: {S == mykeys}; mine minus theirs {len(mykeys - S)}, theirs minus mine {len(S - mykeys)}")
open("out_vmatch.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
