"""u-helpers-ours toy V: an exact small model of the hub's packet-plane helper bound (#942's U8 word: 14 Hamming blocks,
168 face-packet planes of rank 3), in which the least number of helpers is computed exactly and alternatives can be tested.

Model (the hub's own necessary condition, research/structural-gain/factory/README.md, Scholz with Codex, on icekylinx #531):
- the first row whose source support spans packet plane P is created by a two-register addition AT P, so >= 2 registers
  visit P (>= 3 if two retained copies of the packet are needed);
- registers have nested frame chains, and the 168 planes are distinct rank-3 frames, so a register visits at most ONE plane;
- who may visit P: a source s iff its label line u_s lies in P (rule 'src'); a target t, only if output borrowing is
  allowed, iff P lies in t's cap u_t^perp (rule 'cap') or, stricter, iff t is in the hub's compatible-foreign-output list
  for P (rule 'hub', which also asks for the centre frames scattering into t); a helper may visit any plane.
Least helpers = sum_P need_P - (maximum b-matching of sources/targets to planes), computed exactly by max-flow.
This is a LOWER-BOUND model (necessary condition): a count it allows is not shown constructible.
Usage: python3 -I toyV.py PACKET-PLANE-BOUND.json u8-FULL-WORD.json.gz toyV.json"""
import sys, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import engine as E
P = json.load(open(sys.argv[1])); w = E.load(sys.argv[2]); ports = w['ports']; H = 9
planes = [dict(frame=E.basis(pp['frame']), src=pp['source_ids'], hub=pp['compatible_external_outputs']) for pp in P['packet_planes']]
assert len(planes) == 168 and len({p['frame'] for p in planes}) == 168
def inplane(v, f): return E.contained(E.basis([v]), f)
def orth(a, b): return bin(a & b).count('1') % 2 == 0
def maxflow(cap_edges, n, s, t):
    # Edmonds-Karp on an adjacency dict (small graph)
    from collections import deque, defaultdict
    g = defaultdict(dict)
    for u, v, c in cap_edges:
        g[u][v] = g[u].get(v, 0) + c; g[v].setdefault(u, 0)
    flow = 0
    while True:
        par = {s: None}; dq = deque([s])
        while dq and t not in par:
            u = dq.popleft()
            for v, c in g[u].items():
                if c > 0 and v not in par: par[v] = u; dq.append(v)
        if t not in par: return flow
        v = t; bott = 10**9
        while par[v] is not None: bott = min(bott, g[par[v]][v]); v = par[v]
        v = t
        while par[v] is not None: u = par[v]; g[u][v] -= bott; g[v][u] += bott; v = u
        flow += bott
def solve(need, rule_src=True, rule_tgt=None, src_any_plane_with_line=True):
    S, T = 'S', 'T'; edges = []
    for i, p in enumerate(planes): edges.append((('p', i), T, need))
    if rule_src:
        for s in range(E.V):
            edges.append((S, ('s', s), 1))
            for i, p in enumerate(planes):
                if inplane(ports[s], p['frame']): edges.append((('s', s), ('p', i), 1))
    if rule_tgt:
        for t in range(E.V):
            edges.append((S, ('t', t), 1))
            for i, p in enumerate(planes):
                ok = all(orth(ports[t], b) for b in p['frame']) if rule_tgt == 'cap' else (t in p['hub'])
                if ok: edges.append((('t', t), ('p', i), 1))
    f = maxflow(edges, None, S, T); return dict(need_total=need * 168, covered_by_borrowed=f, least_helpers=need * 168 - f)
out = {}
# sanity: each plane holds 4 source lines; each plane lies in exactly how many target caps
out['source_lines_per_plane'] = sorted({sum(inplane(ports[s], p['frame']) for s in range(E.V)) for p in planes})
out['target_caps_per_plane'] = sorted({sum(all(orth(ports[t], b) for b in p['frame']) for t in range(E.V)) for p in planes})
for need, name in ((2, 'one copy'), (3, 'two copies')):
    out[f'{name}: sources only (hub bound, control)'] = solve(need)
    out[f'{name}: sources + targets in whose cap the plane lies'] = solve(need, rule_tgt='cap')
    out[f'{name}: sources + hub-compatible foreign outputs'] = solve(need, rule_tgt='hub')
    out[f'{name}: targets only (no source borrowing)'] = solve(need, rule_src=False, rule_tgt='cap')
json.dump(out, open(sys.argv[3], 'w'), indent=1); print(json.dumps(out, indent=1))
