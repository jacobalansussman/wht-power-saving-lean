"""v-w5-twoshared robust.py: parse w5-twoshared's frame-level outputs (esX, rd2X, rd2Xb, rd2Xc, rd2S) and ask, per core:
is some top excluded by a junk J = y_e q_kj whose main term q_kj has LAYER 3?  Only such junks are immune to the
content-level route the frame tests ignore (a main term of v at layer <= 2 can be cancelled inside t_e's content by
injections u_k' -> t_e / x_j -> u_k' -> t_e), and to the skipped k2 == k1 production source (G2).
A core whose exclusion rests only on junks with lambda(q_kj) <= 2 is reported as NOT ROBUST."""
import re, sys, glob, os
R = os.path.join(os.environ.get('W5_REMOTE_OUT', 'remote-out'), '')   # folder of the track's search outputs (not published)
def parse_new(path):
    cores = []; cur = None
    for line in open(path):
        m = re.match(r"seed (\d+) core (\d+)(?: \((\w+)\))?: (?:sigma .*?, )?tags (\[.*?\]), layer-3 terms of v (\[.*?\])", line)
        if m:
            cur = dict(id=f"{m.group(1)}/{m.group(2)}", lay3=eval(m.group(5)), dead={0: [], 1: []}, topex={0: False, 1: False}, verdict=None)
            cores.append(cur); continue
        m = re.match(r"\s+top (\d) junk (\(\d, \d\)): (\d+) (?:Lemma-D' )?solutions.*alive (\d+)", line)
        if m and cur:
            if int(m.group(4)) == 0: cur['dead'][int(m.group(1))].append(eval(m.group(2)))
            continue
        m = re.match(r"\s+top (\d): (EXCLUDED|passes)", line)
        if m and cur: cur['topex'][int(m.group(1))] = m.group(2) == 'EXCLUDED'; continue
        m = re.match(r"seed (\d+) core (\d+): (EXCLUDED|OPEN)", line)
        if m and cur and cur['id'] == f"{m.group(1)}/{m.group(2)}": cur['verdict'] = m.group(3)
    return cores
def parse_esX(path, seed):
    cores = []; cur = None
    for line in open(path):
        m = re.match(r"core (\d+): tags (\[.*?\]), layer-3 terms of v (\[.*?\])", line)
        if m:
            cur = dict(id=f"{seed}/{m.group(1)}", lay3=eval(m.group(3)), dead={0: [], 1: []}, topex={0: False, 1: False}, verdict=None)
            cores.append(cur); continue
        m = re.match(r"\s+top (\d): (EXCLUDED|passes|OPEN)\s*\[(.*)\]", line)
        if m and cur:
            e = int(m.group(1)); cur['topex'][e] = m.group(2) == 'EXCLUDED'
            for kj, n in re.findall(r"(\(\d, \d\)):(\d+)", m.group(3)):
                if int(n) == 0: cur['dead'][e].append(eval(kj))
            if 'excluded' in line: cur['verdict'] = 'EXCLUDED'
    return cores
groups = {
  'esX (cross 501-504)': [c for s in (501, 502, 503, 504) for c in parse_esX(R + f'ssh/out_esX_{s}.txt', s)],
  'rd2X (cross 501/5, 504/3)': parse_new(R + 'ssh5/out_rd2X_501.txt') + parse_new(R + 'ssh5/out_rd2X_504.txt'),
  'rd2Xb (cross 601-604)': parse_new(R + 'ssh2/out_rd2Xb_A.txt') + parse_new(R + 'ssh2/out_rd2Xb_B.txt'),
  'rd2Xc (cross 31-34 replay)': sum((parse_new(f) for f in [R + 'ssh4/out_rd2Xc_31.txt', R + 'ssh4/out_rd2Xc_32.txt', R + 'ssh5/out_rd2Xc_33.txt', R + 'ssh5/out_rd2Xc_34.txt']), []),
  'rd2S sig1': sum((parse_new(f) for f in [R + 'ssh5/out_rd2S_sig1_a.txt', R + 'ssh5/out_rd2S_sig1_b.txt', R + 'ssh2/out_rd2S_sig1_42.txt', R + 'ssh2/out_rd2S_sig1_44.txt']), []),
  'rd2S sig0': sum((parse_new(f) for f in [R + 'ssh/out_rd2S_sig0_a.txt', R + 'ssh/out_rd2S_sig0_b.txt']), []),
}
for name, cores in groups.items():
    seen = {}
    for c in cores: seen[c['id']] = c          # later duplicates (restarts) override
    rob = []; nonrob = []; notex = []
    for cid, c in seen.items():
        ex = c['topex'][0] or c['topex'][1]
        if not ex: notex.append(cid); continue
        robust = any(c['topex'][e] and any(kj in c['lay3'] for kj in c['dead'][e]) for e in (0, 1))
        (rob if robust else nonrob).append(cid)
    print(f"{name}: {len(seen)} distinct cores; excluded {len(rob) + len(nonrob)} (robust via a layer-3-q junk: {len(rob)}, "
          f"resting only on junks with lambda(q) <= 2: {len(nonrob)} {nonrob}); not excluded here: {len(notex)} {notex}")
