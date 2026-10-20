#!/usr/bin/env python3
"""k-combine (from m-gen's roundtrip.py): parse the GENERATED Lean data module back (x-mutable's lean_raw_to_gcert.py,
text parse only), check it equals the input certificate field by field, then run on the PARSED data the exact mirrors
of every kernel term of the stacked route: scalCheckZ x4, yCheckZ, eCheckZ (zmirror.suite at the preflight's SPar /
E7 SPar) and the 42 zero checks GS.zCheckZ on the y-free twin (zdesignZ.zCheck, rows from the design file).
usage: roundtrip_xz.py <S0.lean> <cert.gcert.json.gz> <scratch.json.gz> <preflight.json> <zdesign.json>"""
import sys, os, json, gzip, contextlib, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm
import lean_raw_to_gcert as lr

def main(s0, cert, tmp, pre, zdp):
    with contextlib.redirect_stdout(sys.stderr):
        lr.main(s0, tmp)
    p = json.loads(gzip.open(tmp, "rt").read())
    c = json.loads(gzip.open(cert, "rt").read())
    keys = ["h", "v", "R", "ports", "frames", "start", "final", "ext", "ret", "scat", "A", "B"]
    diff = [k for k in keys if p[k] != c[k]]
    out = dict(parsed_equals_certificate=not diff, differing=diff)
    P = json.load(open(pre)); sp = tuple(P["SPar"]); ep = tuple(P["epar"])
    lc = zm.lean_data(p); lc["ext"] = p["ext"]
    s = zm.suite(lc, sp, 4, ep, 80, True)
    out["spar"], out["epar"] = sp, ep
    out["kernel_terms_on_parsed_text"] = dict(all_true=s["all_true"], scal=[x["value"] for x in s["scal"]],
                                              ychk=s["ychk"]["value"], echk=[x["value"] for x in s["echk"]])
    # the zero checks on the y-free twin of the PARSED data
    v = p["v"]; isy = lambda r: v <= r < 2 * v
    tw = dict(p); B = []
    for g in p["B"]:
        if g[0] == "out": B.append(["out", g[1], g[2], [x for x in g[3] if not isy(x[0])], g[4]])
        else: B.append(["in", g[1], g[2], [] if isy(g[2]) else g[3]])
    tw["B"] = B; tw["scat"] = {"table": [[] for _ in range(v)]}
    spec = importlib.util.spec_from_file_location("zdz_mirror", os.path.join(os.path.dirname(os.path.abspath(__file__)), "zcheckz_mirror.py"))
    zz = importlib.util.module_from_spec(spec); spec.loader.exec_module(zz)
    zd = json.load(open(zdp)); P5 = tuple(zd["spar"])
    lt = zm.lean_data(tw)
    nr = 2 * v + p["R"]; el = sorted(c["exits"])
    rowsXS = list(range(v)) + list(range(2 * v, nr))
    bitB = lambda l, b, val: [e for e in l if ((l.index(e) >> b) & 1) == val]
    batches = [(list(range(v, v + k)) + bitB(zd["light"], b, val), rowsXS) for b, val, k in zd["bits"]] + \
              [(list(range(v, v + k)) + [e], rowsXS) for e, k in zd["heavy_batches"]] + [(ch, el) for ch in zd["chunks"]]
    vals = [zm.evaluate(zz.zCheck, lt, P5, g, r)["value"] for g, r in batches]
    out["zero_checks_on_parsed_twin"] = dict(n=len(vals), all_true=all(vals))
    print(json.dumps(out))
    assert not diff and s["all_true"] and all(vals)

if __name__ == "__main__":
    main(*sys.argv[1:6])
