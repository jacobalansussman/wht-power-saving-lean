# cert_readable.py CERT.json -- print a certificate of cert_m.py in readable form (pf-duality, 14 Oct 2026)
import sys, json
cert = json.load(open(sys.argv[1])); m = cert["m"]; n = 1 << m
subs = {frozenset({0})}; todo = [frozenset({0})]
while todo:
    nxt = []
    for S in todo:
        for v in range(1, n):
            T = frozenset(S | {x ^ v for x in S})
            if T not in subs: subs.add(T); nxt.append(T)
    todo = nxt
subs = sorted(subs, key=lambda S: (len(S), sorted(S)))
def bits(v): return ''.join(str((v >> b) & 1) for b in range(m))
def nm(i):
    S = subs[i]; d = len(S).bit_length() - 1
    if d == 0: return "L0"
    if d == m: return "M0"
    B = []; span = {0}
    for v in sorted(S):
        if v not in span: B.append(v); span |= {x ^ v for x in span}
    return "N<" + ",".join(bits(v) for v in B) + ">"
def st(fs): return "{" + ", ".join(nm(i) for i in fs) + "}"
for k, r in enumerate(cert["rows"], 1):
    t = r["type"]
    if t == "sub":
        S, a, b = r["witness"]
        print(f"{k:3d}. {r['mult']:>8} x submod  S={st(S)}, a={nm(a)}, b={nm(b)}:  r(S+a)+r(S+b)-r(S+a+b)-r(S) >= 0")
    elif t == "mono":
        S, a = r["witness"]
        print(f"{k:3d}. {r['mult']:>8} x monot   S={st(S)}, a={nm(a)}:  r(S+a)-r(S) >= 0")
    elif t == "ing":
        L = [st(l) for l in r["letters"]]
        print(f"{k:3d}. {r['mult']:>8} x Ingleton A,B,C,D = {', '.join(L)}; T={st(r['T'])}")
    else:
        L = [st(l) for l in r["letters"]]
        print(f"{k:3d}. {r['mult']:>8} x DFZ(65) A,B,C,W,X,Y,Z = {', '.join(L)}; T={st(r['T'])}")
for e in cert["eqs"]:
    d = e["def"]
    if d[0] == "single": print(f"     {e['mult']:>8} x [r({{{nm(d[2])}}}) - r({{{nm(d[1])}}}) = 0]")
    elif d[0] == "pair":
        (s, t), (s0, t0) = d[1], d[2]
        print(f"     {e['mult']:>8} x [delta({nm(s)} over {nm(t)}) - delta({nm(s0)} over {nm(t0)}) = 0]")
    else: print(f"     {e['mult']:>8} x [mean delta_1 = 1]   (the bound)")
