#!/usr/bin/env python3
"""x-mutable driver.  usage: run_checks.py <gx.py> <shea-copies.raw.json.gz> <out dir>
1. control: Shea's complex588 certificate (112 source copies) as rebuilt from his Lean data: histogram, classic gx
   checks, x/y replays, and the classic engine's words (dirt gate + slot-only undo) on random states.
2. the same word WITHOUT the copies (each copy helper q_r = 728 + r folded back into its source x_r, the copy gate
   deleted): classic checker, mutable checker, replays, words with the classic and with the mutable-source undo.
3. counterexample controls: a target read into a helper (breaks hy), a dropped target gate (breaks hid).
"""
import sys, json, gzip, copy, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gxm

def load(p): return json.loads(gzip.open(p, "rt").read())

def uncopy(c):
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

def report(name, c, gx, words=True, trials=2):
    t0 = time.time(); out = {"name": name, "v": c["v"], "R": c["R"]}
    H, per = gxm.hist(c, gx); out["hist"] = H; out["hist_by_class"] = per
    out["N"] = sum(r * n for r, n in H.items())
    out["kinds"] = {ph: dict(k) for ph, k in gxm.kinds(c).items()}
    cb = gxm.structure(c, mutable=False); mb = gxm.structure(c, mutable=True)
    out["classic_structure"] = "pass" if not cb else "FAIL: %d violations, first %s" % (len(cb), cb[0])
    out["mutable_structure"] = "pass" if not mb else "FAIL: %d violations, first %s" % (len(mb), mb[0])
    out["hid"], out["hx_x_columns"] = gxm.hid_hx(c)
    out["hy"] = gxm.hy(c)
    out["hx_all_columns"] = gxm.hx_full(c)
    try:
        cc = dict(c); gx.q2(cc); out["E7_gauged_slots_never_reach_targets"] = True
    except AssertionError as e:
        out["E7_gauged_slots_never_reach_targets"] = "FAIL " + str(e)
    if words:
        out["words_classic_slot_undo"] = gxm.word_tests(c, trials, 7, undo="s")
        out["words_mutable_z_undo"] = gxm.word_tests(c, trials, 7, undo="z")
    out["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(out)); sys.stdout.flush()
    return out

def main(gxp, src, outdir):
    gx = gxm.load_gx(gxp)
    c = load(src)
    res = [report("shea_c588_with_112_source_copies (control)", c, gx)]
    d = uncopy(c)
    gzip.open(os.path.join(outdir, "c588-no-copies.gcert.json.gz"), "wt").write(json.dumps(d))
    res.append(report("c588_no_source_copies (mutable sources)", d, gx))
    # counterexample controls on the no-copy word
    e = copy.deepcopy(d); v = e["v"]
    # turn one s -> y gate of B around (same frame, so the labels stay legal): a helper now READS a target
    for i, g in enumerate(e["B"]):
        if g[0] == "out" and g[2] >= 2 * v and len(g[3]) == 1 and v <= g[3][0][0] < 2 * v and not g[4]:
            t, a, b = g[3][0]; e["B"][i] = ["out", g[1], t, [[g[2], a, b]], []]; break
    res.append(report("control: one s->y gate reversed to y->s (breaks hy)", e, gx, trials=2))
    f = copy.deepcopy(d)
    # drop the first B gate that writes a target from a source (x -> y)
    for i, g in enumerate(f["B"]):
        if g[0] == "out" and g[2] < v and any(v <= t < 2 * v for t, _, _ in g[3]):
            del f["B"][i]; break
    res.append(report("control: one x->y gate deleted (breaks hid)", f, gx, trials=2))
    json.dump(res, open(os.path.join(outdir, "checks.json"), "w"), indent=1)

if __name__ == "__main__":
    main(*sys.argv[1:4])
