#!/usr/bin/env python3
"""g-fields I5, final form: relaxed H2 visit toy restricted to planes that HAVE consumers (both ends of a pair), for
B-phase consumers and for all consumers (compensation/cut reads included). Greedy b-matching (a lower estimate of the
room). Prints the helper visit count left (397 = 392 plane helpers + 5 centre-carrying sources).
Usage: python3 -I relaxH2_consumed.py   (needs timeline.json, consumers.json)"""
import json
from collections import Counter, defaultdict
T=json.load(open('timeline.json')); R=json.load(open('consumers.json'))
ev={int(q):E for q,E in T['ev'].items()}; cls={int(q):k for q,k in T['cls'].items()}; ports=T['ports']; v=112
dot=lambda a,b: bin(a&b).count('1')&1
for label, phase_ok in [('B-phase consumers', lambda t: 11003<=t<12921), ('all consumers', lambda t: True)]:
    cons=defaultdict(set)
    for i,r in enumerate(R):
        for q,k,tv in r['visitors']:
            for (t,ph,role,p,f) in ev[q]:
                if t>tv and role=='r' and cls[p]=='y' and phase_ok(t): cons[i].add(p)
    live=[i for i in range(168) if cons[i]]
    E=[(a,b) for a in live for b in live if a!=b and all(not dot(x,ports[t-v]) for t in cons[b] for x in R[a]['plane'])]
    slots={i:3-sum(1 for q,k,tv in R[i]['visitors'] if k=='x') for i in range(168)}
    und=defaultdict(set)
    for a,b in E: und[a].add(b); und[b].add(a)
    used=0
    for a,b in sorted({tuple(sorted(e)) for e in E}, key=lambda e: len(und[e[0]])+len(und[e[1]])):
        if slots[a]>0 and slots[b]>0: slots[a]-=1; slots[b]-=1; used+=1
    print(label, ': consumed planes', len(live), 'compatible ordered pairs', len(E), 'greedy pair-slots', used, '-> helper visits', 397-used, 'of 397')
