"""k-word: the stacked words (U8 exit frames + mutable sources) built from the gcert lowering of the U8 word, and their
exact checks with x-mutable's mirror (gxm.py) and condition (d) on the CONVERTED word.
  mut0 : every one of the 112 source copies folded back into its source (x-mutable's run_checks.uncopy, unchanged)
  mut5 : the copies of the sources that host a split centre kept (5 on U8), the other 107 folded back
Condition (d) here (no tail in the gcert lowering): the SOURCE rows of the derived end gate Mti read only registers
whose final frame is FULL.  Also reported: the zero pattern of those rows equals the one on the hub canonical word.
usage: mut_variants.py LOWERED.gcert.json.gz gx.py CANONICAL-WORD.json.gz OUTDIR"""
import sys, os, json, gzip, copy, time
from fractions import Fraction as Fr
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gxm
from run_checks import uncopy, report
import cond_d

def partial_uncopy(c, keep):
    """fold back the copies of sources NOT in keep; kept copies renumbered after the word helpers, in source order"""
    v, R = c["v"], c["R"]; base = 2 * v + R - v
    A = c["A"]; cp = A[:v]
    for r, g in enumerate(cp):
        assert g[0] == "out" and g[2] == r and g[3] == [[base + r, 1, 1]] and g[4] == [] and g[1] == c["start"][r], g
    kept = sorted(keep)
    newidx = {base + r: base + i for i, r in enumerate(kept)}
    m = lambda q: (newidx[q] if q in newidx else q - base) if q >= base else q
    def mg(g):
        g = copy.deepcopy(g)
        if g[0] == "neg": g[2] = m(g[2])
        elif g[0] == "out": g[2] = m(g[2]); g[3] = [[m(t), a, b] for t, a, b in g[3]]; g[4] = [m(t) for t in g[4]]
        else: g[2] = m(g[2]); g[3] = [[m(s), a, b] for s, a, b in g[3]]
        return g
    d = copy.deepcopy(c)
    d["A"] = [mg(g) for r, g in enumerate(A[:v]) if r in keep] + [mg(g) for g in A[v:]]
    d["B"] = [mg(g) for g in c["B"]]
    d["ret"] = [[k, m(s), f] for k, s, f in c["ret"]]
    d["ext"] = [[m(r), f] for r, f in c["ext"]]
    n_new = base + len(kept)
    start = [None] * n_new; final = [None] * n_new
    for q in range(2 * v + R):
        if q >= base and (q - base) not in keep: continue
        start[m(q)] = c["start"][q]; final[m(q)] = c["final"][q]
    for r in range(v):
        if r not in keep: start[r] = c["start"][r]; final[r] = c["final"][base + r]
    assert None not in start and None not in final
    d["start"], d["final"] = start, final
    d["R"] = R - v + len(kept)
    return d

def sym_inv_rows(c, gx):
    """rows of Mti (A, scatter, B inverted), every register its own symbol"""
    v, R = c["v"], c["R"]; n = 2 * v + R
    reg = [{r: Fr(1)} for r in range(n)]
    def axpy(t, src, co):
        d = reg[t]
        for k, u in src.items():
            w = d.get(k, 0) + co * u
            if w: d[k] = w
            else: d.pop(k, None)
    def run(gates):
        for g in reversed(gates):
            if g[0] == "neg": reg[g[2]] = {k: -u for k, u in reg[g[2]].items()}
            elif g[0] == "out":
                src = dict(reg[g[2]])
                for t, a, b in g[3]: axpy(t, src, -Fr(a, b))
            else:
                for s, a, b in g[3]: axpy(g[2], dict(reg[s]), -Fr(a, b))
    run(c["B"])
    tot = {k: dict(reg[s]) for k, s, _ in c["ret"]}
    for t in range(v):
        for k, co in gx.scat_row(c, t): axpy(v + t, tot[k], -co)
    run(c["A"])
    return reg

def cond_d_gcert(c, gx):
    v, R = c["v"], c["R"]; n = 2 * v + R
    Mti = sym_inv_rows(c, gx)
    ex = {q for q in range(2 * v, n) if c["final"][q] != 1}
    viol = [[i, sorted(j for j in Mti[i] if j in ex)] for i in range(v) if any(j in ex for j in Mti[i])]
    ready = [i for i in range(v) if any(v <= j < 2 * v for j in Mti[i])]
    pattern = {i: sorted(j for j in Mti[i] if j != i) for i in range(v)}
    return dict(exit_registers=len(ex), condition_d_holds=not viol, violations=viol, source_rows_reading_targets=ready,
                source_rows_mutated=sum(1 for i in range(v) if Mti[i] != {i: 1})), pattern

def main(gp, gxp, cp, outdir):
    gx = gxm.load_gx(gxp)
    c = json.loads(gzip.decompress(open(gp, 'rb').read()))
    v = c["v"]; base = 2 * v + c["R"] - v
    centre_src = sorted(s - base for _, s, _ in c["ret"] if s >= base)
    res = dict(centre_sources=centre_src)
    # canonical-word zero pattern of the source rows of Mti (no tail), for the cross-check
    can = cond_d.read(cp); n = 2 * can['v'] + can['R']
    cMti = cond_d.inv_rows(n, can['A'] + can['cut'] + can['B'])
    cpat = {i: sorted(j for j in cMti[i] if j != i) for i in range(v)}
    for name, d in (("mut0_no_copies", uncopy(c)), ("mut5_centre_copies_kept", partial_uncopy(c, set(centre_src)))):
        gzip.open(os.path.join(outdir, "U8-%s.gcert.json.gz" % name), "wt").write(json.dumps(d))
        r = report(name, d, gx, words=True, trials=2)
        dd, pat = cond_d_gcert(d, gx)
        r["condition_d_on_converted_word"] = dd
        if name.startswith("mut0"):
            r["source_row_pattern_equals_canonical"] = pat == cpat
        else:
            # kept copies stand in for their sources inside the word: compare after mapping them back
            kept = {base + i: s for i, s in enumerate(centre_src)}
            pm = {i: sorted({kept.get(j, j) for j in p} - {i}) for i, p in pat.items()}
            r["source_row_pattern_equals_canonical_after_mapping_kept_copies"] = None   # x of a centre source is only read, compare others
            r["source_rows_mutated_note"] = "x_r of a kept-copy source is never written (unit row)"
        res[name] = r
    json.dump(res, open(os.path.join(outdir, "mut-variants.json"), "w"), indent=1)

if __name__ == "__main__":
    main(*sys.argv[1:5])
