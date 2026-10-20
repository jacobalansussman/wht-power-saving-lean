# decode_dfz.py CERT.json -- name the frames of the DFZ (Fano) rows of a certificate and show the GF(2) picture
# (pf-duality, 14 Oct 2026).  Frame N_s = {(x|0): x in s} + {(0|z): z in s^0} inside F_2^(2m) (PROOF.md Lemma I1).
# For each DFZ row: the letters' images in F_2^(2m) / sum_{t in T} N_t, i.e. what the perfect mod-2 certificate sees.
import sys, json
cert = json.load(open(sys.argv[1])); m = cert["m"]; n = 1 << m
subs = {frozenset({0})}; todo = [frozenset({0})]
while todo:
    nxt = []
    for S in todo:
        for v in range(1, n):
            T = frozenset(S | {x ^ v for x in S})
            if T not in subs: subs.add(T); nxt.append(T)
    todo = nxt
subs = sorted(subs, key=lambda S: (len(S), sorted(S)))
def bits(v): return ''.join(str((v >> b) & 1) for b in range(m))
def basis(S):
    B = []; span = {0}
    for v in sorted(S):
        if v not in span: B.append(v); span |= {x ^ v for x in span}
    return B
def name(i):
    S = subs[i]; d = len(S).bit_length() - 1
    if d == 0: return "L0(s=0)"
    if d == m: return "M0(s=F)"
    return "s=<" + ",".join(bits(v) for v in basis(S)) + ">"
def perp(S): return frozenset(a for a in range(n) if all(bin(a & x).count('1') % 2 == 0 for x in S))
def frame_vecs(i):
    S = subs[i]; return [x for x in S] + [z << m for z in perp(S)]
def rank(vs):
    rows = [];
    for v in vs:
        for r in rows: v = min(v, v ^ r)
        if v: rows.append(v)
    return len(rows)
for row in cert["rows"]:
    if row["type"] != "dfz": continue
    L = [l[0] for l in row["letters"]]; T = row["T"]
    tv = [v for t in T for v in frame_vecs(t)]; rT = rank(tv)
    print(f"DFZ weight {float(eval(row['mult'].replace('/', '/1.0/'))) if False else row['mult']}  T = {[name(t) for t in T]}")
    for lab, f in zip("ABCWXYZ", L):
        print(f"   {lab} = {name(f):24s} raise over T (mod 2) = {rank(tv + frame_vecs(f)) - rT}")
    def r(*fs): return rank(tv + [v for f in fs for v in frame_vecs(f)]) - rT
    A, B, C, W, X, Y, Z = L
    print("   mod-2 dependencies: W in A+B:", r(A, B, W) == r(A, B), " X in A+C:", r(A, C, X) == r(A, C),
          " Y in B+C:", r(B, C, Y) == r(B, C), " Z in A+Y:", r(A, Y, Z) == r(A, Y), " A,B,C independent:", r(A, B, C) == r(A) + r(B) + r(C))
