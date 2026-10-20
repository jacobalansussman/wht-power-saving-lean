"""k-combine (13 Oct 2026): the STACKED certificate U8XM = k-lower's U8X (U8 word of Kothe #869/#870, Chafik #890,
Capalbo #942; exit-frame helpers; the hub tail lowered into phase B; 112 source copies) with the source copies
folded back into their sources (x-mutable's mutable sources), EXCEPT the 5 copies of the sources that host a split
centre (k-word's `partial_uncopy`, copied verbatim; the scatter needs retained totals on slots, gxm rule E2).
Then: block histogram by class, N, cst by an own frame replay (m-gen's prep_cert.blocks, verbatim), the exit-route
invoice identity (E6x), and the exact exit conditions on the total matrix (k-lower's lower_x.total_matrix), now
THREE: exit rows zero off the diagonal, exit columns zero at the other slot rows, exit columns zero at the x rows.
usage: make_cert.py <U8X.gcert.json.gz> <out.gcert.json.gz>"""
import sys, os, json, gzip, copy, hashlib
from collections import Counter
from fractions import Fraction as Q

def partial_uncopy(c, keep):          # k-word/tools/mut_variants.py, verbatim
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

def regs(g):                           # m-gen/tools/prep_cert.py, verbatim
    if g[0] == "out":
        return [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
    return [g[2]] + [s for s, _, _ in g[3]]

def blocks(c):                         # m-gen/tools/prep_cert.py, verbatim
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

def total_matrix(c):                   # k-lower/tools/lower_x.py, verbatim
    v, R = c['v'], c['R']
    n = 2 * v + R
    reg = [{r: Q(1)} for r in range(n)]
    def axpy(t, src, co):
        d = reg[t]
        for k, u in src.items():
            w = d.get(k, 0) + co * u
            if w: d[k] = w
            else: d.pop(k, None)
    def run(gates):
        for g in gates:
            assert g[0] == 'out'
            src = dict(reg[g[2]])
            for t, a, b in g[3]:
                axpy(t, src, Q(a, b))
    run(c['A'])
    tot = {k: dict(reg[s]) for k, s, _ in c['ret']}
    for t in range(v):
        for k, a, b in c['scat']['table'][t]:
            axpy(v + t, tot[k], Q(a, b))
    run(c['B'])
    return reg

def exit_conditions(c, reg):
    v, R = c['v'], c['R']; n = 2 * v + R
    E = set(c['exits'])
    badR = [e for e in sorted(E) if set(reg[e]) - {e}]
    badC = sorted({(i, j) for i in range(2 * v, n) for j in reg[i] if j in E and j != i})
    badX = sorted({(i, j) for i in range(v) for j in reg[i] if j in E})
    yx = all({j: u for j, u in reg[v + t].items() if j < v} == {t: 1} for t in range(v))
    hy = all(not any(v <= j < 2 * v and j != i for j in reg[i]) for i in range(n)) and \
        all(reg[v + t].get(v + t) == 1 for t in range(v))
    xmut = sum(1 for t in range(v) if reg[t] != {t: 1})
    return dict(ExitRows=not badR, ExitRows_bad=badR[:10], ExitCols=not badC, ExitCols_bad=badC[:10],
                ExitColsX=not badX, ExitColsX_bad=badX[:10], hid_y_from_x_identity=yx, hy_y_columns_unit=hy,
                x_rows_not_unit=xmut, exit_diagonals=dict(Counter(str(reg[e].get(e, 0)) for e in E)))

def main(src, dst):
    raw = open(src, 'rb').read()
    c = json.loads(gzip.decompress(raw))
    v, R = c['v'], c['R']; base = 2 * v + R - v
    centre_src = sorted(s - base for _, s, _ in c['ret'] if s >= base)
    d = partial_uncopy(c, set(centre_src))
    assert all(e < base for e in c['exits']) and d['exits'] == c['exits']
    b, N, cst, er = blocks(d)
    h, v, R = d['h'], d['v'], d['R']; n = 2 * v + R
    dim = [len(f) for f in d['frames']]
    inv = [N + er, sum(dim[d['final'][r]] for r in range(2 * v, n)) + 2 * v * (h - 1) + cst]
    assert inv[0] == inv[1], ('E6x', inv)
    d['blocks'], d['N'], d['cst'] = b, N, cst
    d['derived_from'] = ("k-lower's U8X certificate (sha256 %s; U8 word of Kothe #869/#870, Chafik #890, Capalbo #942, "
                         "hub canonical sha256 24db880b...; exit helpers, hub tail in phase B, 112 source copies) with "
                         "the copies folded back except the 5 centre-source copies %s (k-combine 13 Oct 2026)"
                         % (hashlib.sha256(raw).hexdigest(), centre_src))
    out = json.dumps(d, separators=(',', ':'), sort_keys=True).encode()
    open(dst, 'wb').write(gzip.compress(out, mtime=0))
    hist = Counter()
    for w in b.values():
        for r, k in w.items(): hist[int(r)] += k
    ec = exit_conditions(d, total_matrix(d))
    print(json.dumps(dict(src_sha256=hashlib.sha256(raw).hexdigest(), out_sha256=hashlib.sha256(open(dst, 'rb').read()).hexdigest(),
                          json_sha256=hashlib.sha256(out).hexdigest(), centre_sources_kept=centre_src, v=v, R=R, N=N,
                          cst=cst, entrance_rank_sum=er, invoice_E6x=inv, exits=len(d['exits']), blocks=b,
                          hist=dict(sorted(hist.items())), gates=dict(A=len(d['A']), B=len(d['B'])), exit_conditions=ec), indent=1))
    assert ec['ExitRows'] and ec['ExitCols'] and ec['ExitColsX'] and ec['hid_y_from_x_identity'] and ec['hy_y_columns_unit']

if __name__ == '__main__':
    main(*sys.argv[1:3])
