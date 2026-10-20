"""d-mcopy-code STEP 2c: can any EXISTING register (fixed chain; helper, source, or another target) carry any needed
source of an M group to one of the group's targets?  Necessary: frames G1 <= G2 of r's chain with (a) some giver
(helper/source, or r itself; a helper/source r may not read a target) holding needed content at G1, and (b) G2 a
frame of a group target's chain (a common gate frame; for r a target this is a target-target read).
Timing ignored (relaxation).  Also checks the packet-frame visitors (only creators + M copy).
Usage: python3 -I existing_relays.py FULL-WORD.json.gz"""
import sys, json
exec(open('chain_cover.py').read().split("gi = int(sys.argv[2])")[0])
Gs = json.load(open('mc_groups.json')); allM = set(q for g in Gs for q in g['members'])
AG = json.load(open('all_groups.json'))
start, end, chain, paid, runs = frames(w)
prog = program(w)
T_CLEAN = min(t for t, (ph, i, e) in enumerate(prog) if ph == 'cleanup')
tagtime = {(ph, i): t for t, (ph, i, e) in enumerate(prog)}
def ttime(tag):
    if tag is None: return -1
    if tag[0] == 'final': return 10**9
    return tagtime[tuple(tag)]
total = 0; report = []
for gi, A in enumerate(AG):
    T = A['targets']; mem = set(A['members'])
    rows = [{q: F(1)} for q in range(n)]; need = set()
    for t, (ph, i, e) in enumerate(prog):
        if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
            xs = [s for s in rows[e[2]] if s < v]
            if len(xs) == 8: need |= set(xs)
        apply(e, rows)
    # giver frames (helpers/sources, not M, not targets) holding needed content, windows capped at cleanup
    wins = defaultdict(list)
    for g in range(n):
        if g in allM or v <= g < 2*v: continue
        ch = chain[g]
        for k, (tag, U) in enumerate(ch):
            a = ttime(tag); b = min(ttime(ch[k+1][0]) if k + 1 < len(ch) else 10**9, T_CLEAN)
            wins[g].append((a, b, tuple(U)))
    rows = [{q: F(1)} for q in range(n)]; gf = set()
    def upd(q, t):
        for (a, b, U) in wins.get(q, []):
            if a <= t < b and any(s in need for s in rows[q]): gf.add(U)
    for q in wins: upd(q, -1)
    for t, (ph, i, e) in enumerate(prog):
        apply(e, rows)
        for q in (e[1:3] if e[0] == 'add' else e[1:2]): upd(q, t)
    tframes = set(tuple(U) for t in T for tag, U in chain[t] if U)
    hits = []
    for r in range(n):
        if r in allM or r in T: continue
        ch = [tuple(U) for tag, U in chain[r]]
        for i1, G1 in enumerate(ch):
            if G1 not in gf: continue
            for G2 in ch[i1:]:
                if G2 in tframes: hits.append((r, len(G1), len(G2)))
    total += len(hits)
    report.append(dict(group=gi, giver_frames=len(gf), existing_relay_paths=hits[:20]))
    print('group', gi, 'giver frames with needed content', len(gf), ' existing-register relay paths:', len(hits), hits[:6])
print('TOTAL existing-register relay paths over 14 groups:', total)
json.dump(report, open('existing_relays.json', 'w'), indent=0)
