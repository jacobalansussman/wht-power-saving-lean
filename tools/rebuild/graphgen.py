"""graphgen.py: generator of the paired-cube signed addition circuit.

ADAPTED from two scripts of the outside repository (github.com/CrocSwap/integer-mult-bounds, Apache-2.0) at
commit 4a3c769.  They were read as text and never run.
  class Gen  follows class Graph of scripts/paired_cube/graph.py method by method (add; sum, here total; module,
             here run_module; configured_local_channels, here local_configured; finish), with other names for
             its fields.  Notice of that file: "Copyright 2026 icekylinx. Apache-2.0." and "Developed with GPT-6
             Astra assistance; integrated with Codex assistance."  Its method configured_local_channels came
             with pull request #168 (eumemic).
  Gen.fuse   follows the branch 'f8:' of merge_outputs() of scripts/paired_cube/modules.py.  That file has the
             same notice and marks this part "Output merging (pmerge.py variants; eumemic, Claude assistance;
             Apache-2.0)".
Changed here: every node carries a tag (the dictionary kindof) that says which part of the circuit made it; the
default local circuit and the centre routing by a tree, which this word does not use, are left out; graph()
returns fewer fields.  NOTICE, section 6, lists all files of this kind.
Mathematics: the note notes/paired-cube-construction.tex of the outside repository.  The node numbering is that
of the outside producer, because the frozen data files of the outside repository (matching arcs, moved frames,
reuse pairs) name nodes, operations and roles by index.  Nodes are 0-based here: 0..v-1 are the ports.

ports: p coordinates pairs, h = 2p; a port T = (I, bits), I a 3-subset of [p], vector e_{2i+b_i} summed: weight 3.
value of node x = value(a) + sign * value(b) for args[x] = [a, b].
"""
from itertools import combinations, product

BITS = list(product(range(2), repeat=3))


class Gen:
    def __init__(self, p):
        self.p, self.h = p, 2 * p
        self.cubes = list(combinations(range(p), 3))
        self.labels = [tuple(2 * i + b for i, b in zip(I, bits)) for I in self.cubes for bits in BITS]
        self.v = len(self.labels)
        self.args = [None] * self.v
        self.signs = [1] * self.v
        self.supp = [1 << i for i in range(self.v)]
        self.memo = {}
        self.port = {(I, bits): 8 * j + k for j, I in enumerate(self.cubes) for k, bits in enumerate(BITS)}
        self.F, self.A, self.G = {}, {}, {}
        self.roots, self.centers = [], []
        self.kindof = {}          # node -> role kind tag of the sub-circuit that created it

    def add(self, a, b, sign=1, tag="?"):
        if a is None:
            assert sign == 1
            return b
        if b is None:
            return a
        if sign == 1 and a > b:
            a, b = b, a
        key = (a, b, sign)
        x = self.memo.get(key)
        if x is None:
            assert not self.supp[a] & self.supp[b], key
            x = len(self.args)
            self.memo[key] = x
            self.args.append([a, b])
            self.signs.append(sign)
            self.supp.append(self.supp[a] | self.supp[b])
            self.kindof[x] = tag
        return x

    def total(self, xs, tag):
        xs = list(xs)
        if not xs:
            return None
        while len(xs) > 1:
            xs = [self.add(xs[i], xs[i + 1], 1, tag) if i + 1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
        return xs[0]

    def run_module(self, mod, inputs, tag):
        args, n = mod["args"], mod.get("input_count", len(inputs))
        one = args[0] == [0, 0] and len(args) > n and args[n] == [0, 0]
        img = ([None] if one else []) + list(inputs)
        for a, b in args[(n + 1 if one else n):]:
            img.append(self.add(img[a], img[b], 1, tag))
        return [img[x] for x in mod["roots"]]

    # ---- the 13 local outputs of a cube: A[i,a] (6 face sums), G[j,k,mode] (6 signed edge channels), F (total)
    def local_configured(self, cfg):
        Ac = {tuple(int(c) for c in k.split(",")): w for k, w in cfg["A"].items()}
        Gc = {tuple(int(c) for c in k.split(",")): w for k, w in cfg["G"].items()}
        for I in self.cubes:
            def at(d):
                b = [0] * 3
                for q, w in d.items():
                    b[q] = w
                return self.port[I, tuple(b)]

            def edge(d, fixed):
                lo = at({**fixed, d: 0})
                hi = at({**fixed, d: 1})
                return self.add(lo, hi, 1, "local")

            def signed(xa, xb):
                return (self.add(xa, xb, -1, "local"), 1) if xa < xb else (self.add(xb, xa, -1, "local"), -1)
            for j, k in combinations(range(3), 2):
                r = 3 - j - k
                for mode in range(2):
                    u, w = 0, mode
                    kd = Gc[j, k, mode]
                    if kd == "e":
                        n1 = edge(r, {j: u, k: w})
                        n2 = edge(r, {j: 1 - u, k: 1 - w})
                        node = self.add(n1, n2, -1, "local")
                    else:
                        far = (1, 0) if kd == "l" else (0, 1)
                        t1 = signed(at({j: u, k: w, r: 0}), at({j: 1 - u, k: 1 - w, r: far[0]}))
                        t2 = signed(at({j: u, k: w, r: 1}), at({j: 1 - u, k: 1 - w, r: far[1]}))
                        if t1[1] == 1:
                            node = self.add(t1[0], t2[0], t2[1], "local")
                        else:
                            assert t2[1] == 1
                            node = self.add(t2[0], t1[0], -1, "local")
                    self.G[I, I[j], I[k], mode] = node
            for i in range(3):
                j, k = [q for q in range(3) if q != i]
                for a in range(2):
                    kd = Ac[i, a]
                    if kd == "fd":
                        n1 = self.add(at({i: a, j: 0, k: 0}), at({i: a, j: 1, k: 1}), 1, "local")
                        n2 = self.add(at({i: a, j: 0, k: 1}), at({i: a, j: 1, k: 0}), 1, "local")
                    else:
                        d = int(kd[1])
                        assert kd[0] == "e" and d != i
                        o = j if d == k else k
                        n1 = edge(d, {i: a, o: 0})
                        n2 = edge(d, {i: a, o: 1})
                    self.A[I, I[i], a] = self.add(n1, n2, 1, "local")
            f = cfg["F"]
            self.F[I] = self.add(self.A[I, I[f], 0], self.A[I, I[f], 1], 1, "local")

    def finish(self, triple, pair, allbut, local):
        p = self.p
        self.local_configured(local)
        n0 = len(self.args)
        D = dict(zip(self.cubes, self.run_module(triple, [self.F[I] for I in self.cubes], "triple")))
        P, Q = {}, {}
        for i in range(p):
            others = [a for a in range(p) if a != i]
            prs = list(combinations(others, 2))
            for bit in range(2):
                ins = [self.A[tuple(sorted((i,) + K)), i, bit] for K in prs]
                for K, node in zip(prs, self.run_module(pair, ins, "pair")):
                    P[tuple(sorted((i,) + K)), i, bit] = node
        for i, j in combinations(range(p), 2):
            others = [a for a in range(p) if a not in (i, j)]
            for mode in range(2):
                ins = [self.G[tuple(sorted((i, j, k))), i, j, mode] for k in others]
                for k, node in zip(others, self.run_module(allbut, ins, "allbut")):
                    Q[tuple(sorted((i, j, k))), i, j, mode] = node
        n1 = len(self.args)

        def emit(I, bl, node, num, ch):
            self.roots.append(dict(node=node, targets=[self.port[I, b] for b in bl],
                                   coefficients=["%d/2" % num] * len(bl), kind="side", channel=ch))
        for I in self.cubes:
            emit(I, BITS, D[I], 1, "disjoint")
            for a in range(2):
                emit(I, [b for b in BITS if b[0] == a], P[I, I[0], 1 - a], 1, "face0")
            for a, b in product(range(2), repeat=2):
                emit(I, [(a, b, c) for c in range(2)], P[I, I[1], 1 - b], 1, "face1")
            for a, b in product(range(2), repeat=2):
                emit(I, [(a, b, c) for c in range(2)], Q[I, I[0], I[1], a ^ b], 2 * a - 1, "edge01")
            for bits in BITS:
                emit(I, [bits], P[I, I[2], 1 - bits[2]], 1, "face2")
            for ii, jj in ((0, 2), (1, 2)):
                for bits in BITS:
                    emit(I, [bits], Q[I, I[ii], I[jj], bits[ii] ^ bits[jj]], 2 * bits[ii] - 1, "edge%d%d" % (ii, jj))
        for i in range(p):                       # stars: S_c = sum of the ports through coordinate c = 2i+e
            for e in range(2):
                self.centers.append(self.total((self.A[I, i, e] for I in self.cubes if i in I), "star"))
        for c, node in enumerate(self.centers):
            self.roots.append(dict(node=node, targets=list(range(self.v)),
                                   coefficients=["1/3" if c in t else "-1/6" for t in self.labels],
                                   kind="center", coordinate=c))
        self.counts = dict(local=n0 - self.v, query=n1 - n0, center=len(self.args) - n1)
        return self

    def fuse(self, code):
        """Follows pull request #168 (eumemic), which says of its fusion "This follows #163's same-singleton
        fusion" and credits it so: "Same-singleton output fusion: @chafreaky (#163); the fold orders here come
        from our own search."  Each port's three single-target channels (face2, edge02, edge12) become one signed
        sum, folded in the order perms[code[t % 8]] (merge_outputs of the outside scripts/paired_cube/modules.py,
        variant 'f8:<code>', read as text)."""
        from fractions import Fraction
        perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        names = ("face2", "edge02", "edge12")
        side = [r for r in self.roots if r["kind"] == "side"]
        centre = [r for r in self.roots if r["kind"] != "side"]
        single, keep = {}, []
        for r in side:
            if r["channel"] in names:
                single[r["targets"][0], r["channel"]] = (r["node"], Fraction(r["coefficients"][0]))
            else:
                keep.append(r)
        new = []
        for t in range(self.v):
            parts = [single[t, ch] for ch in names]
            order = perms[int(code[t % 8])]
            node, c0 = parts[order[0]]
            for k in order[1:]:
                nd, c = parts[k]
                assert c / c0 in (1, -1)
                node = self.add(node, nd, int(c / c0), "fuse")
            new.append(dict(node=node, targets=[t], coefficients=[str(c0)], kind="side", channel="fused"))
        self.roots = keep + new + centre
        return self

    def graph(self):
        return dict(p=self.p, h=self.h, v=self.v, labels=self.labels,
                    inputs=[sum(1 << x for x in t) for t in self.labels],
                    args=self.args, signs=self.signs, roots=self.roots, centers=self.centers)
