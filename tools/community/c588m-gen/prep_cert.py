#!/usr/bin/env python3
"""m-gen (13 Oct 2026): the input certificate of the copy-free #588 build, with its provenance checked.

1. Load Ryan Shea's published candidate (lean/complex588/artifacts/candidate.gcert.json.gz of
   shea256/fourier-transform-below-nlogn @ 5239484: Jeremy Kothe's #588 atm3 circuit with 112 source copies).
2. Fold the 112 source copies back into their sources (x-mutable's `uncopy`, copied verbatim from
   x-mutable/tools/run_checks.py): the copy gate x_r -> q_r is deleted and every later use of q_r = 728 + r
   becomes a use of x_r.
3. Check that the result equals x-mutable's out/c588-no-copies.gcert.json.gz field by field (that file was
   rebuilt from Shea's GENERATED Lean data, so this ties both routes together).
4. Recompute the block histogram by register class, N and cst by an own frame replay (no trust in the stored
   fields), compare with the stored ones, and write the certificate with deterministic gzip (mtime 0).

usage: prep_cert.py <shea candidate.gcert.json.gz> <x-mutable c588-no-copies.gcert.json.gz> <out.gcert.json.gz>
Pure data handling: nothing of the inputs is executed.
"""
import sys, json, gzip, copy, hashlib
from collections import Counter


def load(p):
    return json.loads(gzip.open(p, "rt").read())


def uncopy(c):  # x-mutable/tools/run_checks.py, verbatim
    v, R = c["v"], c["R"]; base = 2 * v + R - v      # copies are the LAST v registers
    A = c["A"]
    cp = A[:v]
    for r, g in enumerate(cp):
        assert g[0] == "out" and g[2] == r and g[3] == [[base + r, 1, 1]] and g[4] == [] and g[1] == c["start"][r], g
        assert c["start"][base + r] == 0
    m = lambda q: q - base if q >= base else q
    def mg(g):
        g = copy.deepcopy(g)
        if g[0] == "neg": g[2] = m(g[2])
        elif g[0] == "out": g[2] = m(g[2]); g[3] = [[m(t), a, b] for t, a, b in g[3]]; g[4] = [m(t) for t in g[4]]
        else: g[2] = m(g[2]); g[3] = [[m(s), a, b] for s, a, b in g[3]]
        return g
    d = copy.deepcopy(c)
    d["A"] = [mg(g) for g in A[v:]]
    d["B"] = [mg(g) for g in c["B"]]
    d["ret"] = [[k, m(s), f] for k, s, f in c["ret"]]
    d["ext"] = [[m(r), f] for r, f in c["ext"]]
    assert all(r < base for r, _ in d["ext"])
    d["start"] = c["start"][:base]
    d["final"] = [c["final"][base + r] if r < v else c["final"][r] for r in range(base)]
    d["R"] = R - v
    return d


def regs(g):
    if g[0] == "out":
        return [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    return [g[2]] + [s for s, _, _ in g[3]]


def blocks(c):
    """block histogram by class (x, y, s, c = centre copies at the scatter), as gxcore.chk_replay counts it"""
    v, R = c["v"], c["R"]
    dim = [len(f) for f in c["frames"]]
    cls = lambda r: "x" if r < v else "y" if r < 2 * v else "s"
    H = {k: Counter() for k in "xysc"}
    cur = list(c["start"])

    def climb(r, f):
        if cur[r] != f:
            H[cls(r)][dim[f] - dim[cur[r]]] += 1
            cur[r] = f
    for g in c["A"]:
        for r in regs(g):
            climb(r, g[1])
    for _, s, f in sorted(c["ret"]):
        assert cur[s] == f, "retained slot not at its frame"
        H["c"][dim[f]] += 1
    for g in c["B"]:
        for r in regs(g):
            climb(r, g[1])
    for r in range(2 * v + R):
        climb(r, c["final"][r])
    out = {k: {str(r): n for r, n in sorted(w.items())} for k, w in H.items()}
    N = sum(r * n for w in H.values() for r, n in w.items())
    cst = sum(r * n for r, n in H["c"].items())
    ext_rank = sum(dim[f] for _, f in c["ext"])
    return out, N, cst, ext_rank


def main(src, xm, dst):
    raw = open(src, "rb").read()
    c = json.loads(gzip.decompress(raw))
    b0, N0, cst0, e0 = blocks(c)
    assert b0 == c["blocks"] and N0 == c["N"] and cst0 == c["cst"], "Shea's stored histogram differs from the replay"
    d = uncopy(c)
    x = load(xm)
    keys = ["format", "h", "v", "R", "ports", "frames", "start", "final", "ext", "ret", "scat", "A", "B"]
    diff = [k for k in keys if d[k] != x[k]]
    assert not diff, "differs from x-mutable's certificate in %s" % diff
    b, N, cst, er = blocks(d)
    h, v, R = d["h"], d["v"], d["R"]
    assert N + er == R * h + 2 * v * (h - 1) + cst, "E6 identity"
    for k in ("blocks", "N", "cst"):
        assert x[k] == {"blocks": b, "N": N, "cst": cst}[k], "x-mutable's stored %s differs from the replay" % k
    d["blocks"], d["N"], d["cst"] = b, N, cst
    d["derived_from"] = ("Ryan Shea's complex588 candidate (shea256/fourier-transform-below-nlogn @ 5239484, "
                         "lean/complex588/artifacts/candidate.gcert.json.gz, sha256 %s: Jeremy Kothe's #588 atm3 with "
                         "112 source copies) with the copies folded back into the sources (x-mutable's uncopy; "
                         "m-gen 13 Oct 2026)" % hashlib.sha256(raw).hexdigest())
    out = json.dumps(d, separators=(",", ":"), sort_keys=True).encode()
    open(dst, "wb").write(gzip.compress(out, mtime=0))
    hist = Counter()
    for w in b.values():
        for r, n in w.items():
            hist[int(r)] += n
    print(json.dumps(dict(shea_candidate_sha256=hashlib.sha256(raw).hexdigest(), shea=dict(R=c["R"], N=N0, blocks=b0),
                          equal_to_x_mutable=True, out=dst, out_sha256=hashlib.sha256(open(dst, "rb").read()).hexdigest(),
                          json_sha256=hashlib.sha256(out).hexdigest(), v=v, R=R, N=N, cst=cst, entrance_rank_sum=er,
                          blocks=b, hist=dict(sorted(hist.items())), gates=dict(A=len(d["A"]), B=len(d["B"])))))


if __name__ == "__main__":
    main(*sys.argv[1:4])
