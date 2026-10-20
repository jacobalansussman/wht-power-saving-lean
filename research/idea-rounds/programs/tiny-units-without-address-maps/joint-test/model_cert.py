#!/usr/bin/python3
"""joint/cert/model_cert.py -- agent "cert", 10 Oct 2026.  My own computation of capacities and demands at block size 3,
f = 2, from 8 x 8 complex matrices built here (pure Python, exact: all entries are dyadic).  No import from any other
agent's program.  Writes model_cert.txt (7 x 15 x 15 capacities, then 15 x 15 demands) for certfind.c and certcheck.c,
and, given a second file, compares with another agent's model (model_m3.json, not published).
My block numbering: 0 = (x = 0); x = 1..7 -> block x = (x, bit 0); block 7 + x = (x, bit 1); bit = x.z."""
import sys, json
n = 8
CA = (1 + 1j) / 2; CB = (1 - 1j) / 2
def dot(a, b): return bin(a & b).count("1") & 1
def mm(a, b): return [[sum(a[i][k] * b[k][j] for k in range(n) if a[i][k] != 0) for j in range(n)] for i in range(n)]
def tr(a): return [[a[j][i] for j in range(n)] for i in range(n)]
def shift(z): return [[1 + 0j if j == (i ^ z) else 0j for j in range(n)] for i in range(n)]
def dirm(z):
    s = shift(z); return [[CA * (i == j) + CB * s[i][j] for j in range(n)] for i in range(n)]
def phase(z, c): return [[(1j if ((dot(z, i) + c) & 1) else 1 + 0j) if i == j else 0j for j in range(n)] for i in range(n)]
def pauli(x, z):
    s = shift(x); return [[s[i][j] * (-1.0) ** dot(z, j) for j in range(n)] for i in range(n)]
KER = mm(mm(dirm(1), dirm(2)), dirm(4))
C2 = [[CA, CB], [CB, CA]]
assert all(KER[i][j] == C2[(i >> 2) & 1][(j >> 2) & 1] * C2[(i >> 1) & 1][(j >> 1) & 1] * C2[i & 1][j & 1] for i in range(n) for j in range(n)), "kernel is not C x C x C"
labs = [(x, z) for x in range(n) for z in range(n)]
idx = {l: i for i, l in enumerate(labs)}
PA = [pauli(x, z) for x, z in labs]
def perm_of(g):
    """M -> g M g^T on Pauli matrices: returns the list p with g P_q g^T = (power of i) * P_p[q]; asserts it is so."""
    gt = tr(g); out = []
    for q in range(64):
        img = mm(mm(g, PA[q]), gt)
        nz = [j for j in range(n) if img[0][j] != 0]
        assert len(nz) == 1; x = nz[0]
        c = img[x][0]; assert c in (1, 1j, -1, -1j), c
        z = 0
        for kb in range(3):
            e = 1 << kb; r = img[x ^ e][e] / c
            assert r in (1, -1)
            if r == -1: z |= e
        P = PA[idx[(x, z)]]
        assert all(img[i][j] == c * P[i][j] for i in range(n) for j in range(n)), "not a Pauli matrix times a power of i"
        out.append(idx[(x, z)])
    assert sorted(out) == list(range(64))
    return out
free = [perm_of(shift(z)) for z in range(1, n)] + [perm_of(phase(z, c)) for z in range(n) for c in (0, 1)]
DIR = {w: perm_of(dirm(w)) for w in range(1, n)}
KP = perm_of(KER)
par = list(range(64))
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for p in free:
    for q in range(64): par[find(q)] = find(p[q])
roots = sorted(set(find(q) for q in range(64)))
assert len(roots) == 15
def bid(q):
    x, z = labs[q]; return 0 if x == 0 else (x if dot(x, z) == 0 else 7 + x)
for r in roots:
    assert len(set(bid(q) for q in range(64) if find(q) == r)) == 1
assert len(set(bid(r) for r in roots)) == 15
blk = [bid(q) for q in range(64)]
size = [blk.count(b) for b in range(15)]
cap = {w: [[0] * 15 for _ in range(15)] for w in DIR}
for w, p in DIR.items():
    for q in range(64): cap[w][blk[q]][blk[p[q]]] += 1
dem = [[0] * 15 for _ in range(15)]
for q in range(64): dem[blk[q]][blk[KP[q]]] += 1
# closed forms, checked on all 64 labels
assert all(DIR[w][idx[(x, z)]] == idx[((x ^ w) if dot(w, z) == 0 else x, z)] for w in DIR for x, z in labs)
assert all(KP[idx[(x, z)]] == idx[(x ^ 7 ^ z, z)] for x, z in labs)
with open(sys.argv[1], "w") as fh:
    for w in range(1, n):
        for s in range(15): fh.write(" ".join(map(str, cap[w][s])) + "\n")
    for s in range(15): fh.write(" ".join(map(str, dem[s])) + "\n")
print("15 blocks, sizes", size)
print("labels a block receives from other blocks under the kernel:", [sum(dem[s][d] for s in range(15) if s != d) for d in range(15)])
print("closed forms hold on all 64 labels: dir(w): (x|z) -> (x+w|z) if w.z = 0 else fixed; kernel: (x|z) -> (x + not z | z)")
if len(sys.argv) > 2:
    M = json.load(open(sys.argv[2])); name = [tuple(v) for v in M["name"]]
    mine = lambda t: 0 if t[0] == 0 else (t[0] if t[1] == 0 else 7 + t[0])
    bc = sum(cap[w][mine(s)][mine(d)] != M["cap"][str(w)][i][j] for w in range(1, n) for i, s in enumerate(name) for j, d in enumerate(name))
    bd = sum(dem[mine(s)][mine(d)] != M["dem"][i][j] for i, s in enumerate(name) for j, d in enumerate(name))
    print("compared with model_m3.json: %d capacity entries, %d differ; %d demand entries, %d differ; sizes equal: %s" % (7 * 225, bc, 225, bd, [size[mine(s)] for s in name] == M["size"]))
