"""d-mcopy-code variant: existing registers that could become the third visitor of a packet frame F_j by splitting ONE
run (inserting F_j between consecutive chain frames c_i < F_j < c_{i+1}) and that later reach U_t of a reader t of P_j
(the only frame of t's chain containing F_j).  Timing: r must be at c_i at some time while the holder is at F_j.
Usage: python3 -I split_candidates.py FULL-WORD.json.gz"""
import sys, json
exec(open('chain_cover.py').read().split("gi = int(sys.argv[2])")[0])
Gs = json.load(open('mc_groups.json')); AG = json.load(open('all_groups.json'))
allM = set(q for g in Gs for q in g['members'])
start, end, chain, paid, runs = frames(w)
prog = program(w)
tagtime = {(ph, i): t for t, (ph, i, e) in enumerate(prog)}
def ttime(tag):
    if tag is None: return -1
    if tag[0] == 'final': return 10**9
    return tagtime[tuple(tag)]
tot = Counter(); ex = []
for gi, A in enumerate(AG):
    for p in A['packets']:
        if p['nsrc'] != 8: continue
        Fj = tuple(p['frame']); hold = p['holder']
        hc = chain[hold]; k = [tuple(U) for tag, U in hc].index(Fj)
        ha, hb = ttime(hc[k][0]), (ttime(hc[k+1][0]) if k + 1 < len(hc) else 10**9)
        for t in p['readers']:
            Ut = tuple(end[t])
            for r in range(n):
                if r in allM or r == t: continue
                c = chain[r]; fr_ = [tuple(U) for tag, U in c]
                if Ut not in fr_: continue
                for i in range(len(c) - 1):
                    a, b = c[i][1], c[i+1][1]
                    if contained(a, Fj) and contained(Fj, b) and tuple(a) != Fj and tuple(b) != Fj:
                        if fr_.index(Ut) < i + 1: continue
                        ra, rb = ttime(c[i][0]), ttime(c[i+1][0])
                        overlap = max(ra, ha) < min(rb, hb)
                        cl = 'x' if r < v else ('y' if r < 2*v else 'h')
                        tot[(cl, overlap)] += 1
                        if overlap and len(ex) < 12: ex.append(dict(group=gi, packet_M=p['M'], reader=t, r=r, cls=cl, run=(len(a), len(b)), r_window=(ra, rb), holder_window=(ha, hb)))
print('split-run third-visitor candidates (class, time-overlap with holder at F_j):', dict(tot))
for e in ex: print('  ', e)
json.dump(dict(counts={f'{a}:{b}': x for (a, b), x in tot.items()}, examples=ex), open('split_candidates.json', 'w'), indent=0)
