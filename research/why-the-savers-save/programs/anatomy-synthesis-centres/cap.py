"""an-synth calc/cap.py -- one formula for all three savers, the design x template cap table, the
two-sided-slack prediction, and what the 1e-1 bar would need. Our own code, stdlib, < 1 s.
Run: python3 -I calc/cap.py        (reads calc/ade.json written by calc/ade.py)

First-order law (an-hub's a ~ D/E, an-ours' D/N, an-openai's identity M = Wm - 2N + 2L):
    a  ~  D / E,   D = mW - (paid rank),   E = sum over paid children r ln(m/r)
       =  (D/(W m)) / lbar,   lbar = E / (paid rank) ~ mean of ln(m/r)
with D/W = rho * (u - c/w): rho = data families / all arrays, u = free coordinates per family,
c/w = centre (bridge) rank per family.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))

def first_order(D, W, m, lbar):
    return D / (W * m * lbar)

print("== (a) calibration: one first-order law against each source's own figure ==")
rows = [
    # name, D, W, m, lbar, reported a, source of the inputs
    ("OpenAI explicit DFT unit (h=100), per-dimension pricing",
     6871402692000000, 1873807244643542670000, 10**6, math.log(10**6), 2.65e-13, "an-openai calc/family.py"),
    ("hub T74 (#1839 line), per call: w=99.556, D=35.556, E: r_eff 3.945",
     35.5556, 99.5556, 73, math.log(73 / 3.9448), 1.6796718e-3, "an-hub out/anatomy.txt"),
    ("our U8X0HQ (Lean-checked), D=2560, W=7491, lbar 2.885",
     2560, 7491, 73, 2.885, 1.6256629e-3, "an-ours run1.out"),
]
for name, D, W, m, lb, rep, src in rows:
    a = first_order(D, W, m, lb)
    print(f"  {name:70s} D/(Wm)={D/(W*m):.4e} lbar={lb:.3f}  a~{a:.4e}  reported {rep:.4e}  ({src})")

print()
print("== (b) design x template: D/(Wm) of the hub's template m = 1 + (k-1)(2P+1), helpers ignored ==")
print("   per run: D = 2Pv - (2P+1) cst ; one free coordinate per family (u = 1); cst = t*(k-1) on v = n - t ports")
print("   (the frame-line-removed variant that gives P112: v = 112, cst = 64) ; best P over 1..8")
ade = json.load(open(os.path.join(HERE, "ade.json")))
designs = []
for d in ade:
    if d["design"].startswith("E8 as 60"):
        continue
    designs.append((d["design"], d["n"], d["t"], d["window_k"], d["junk_per_line"]))
for h in (9, 12, 16, 20, 25, 30):  # Alon triples: n = C(h,3) ports, t = h+1 totals, window h
    n = math.comb(h, 3)
    s = math.comb(h - 3, 3) + 3 * (h - 3)
    designs.append((f"Alon triples h={h} (OpenAI's design)", n + h + 1, h + 1, h + 1, s))

def template(n, t, k):
    v = n - t
    cst = t * (k - 1)
    best = None
    for P in range(1, 9):
        m = 1 + (k - 1) * (2 * P + 1)
        D = 2 * P * v - (2 * P + 1) * cst
        x = D / (2 * P * v * m)
        if best is None or x > best[0]:
            best = (x, P, m, D / (2 * P * v))
    first_m = None
    for P in range(1, 40):
        if 2 * P * v - (2 * P + 1) * cst > 0:
            first_m = 1 + (k - 1) * (2 * P + 1); break
    return cst / v, best, first_m

print(f"   {'design':40s} {'ports v':>7s} {'cst/v':>6s} {'junk/port':>9s} {'best P':>6s} {'m':>4s} {'D/W':>6s} {'D/(Wm)':>9s} {'first m saving':>14s}")
table = []
for name, n, t, k, s in designs:
    ratio, best, first_m = template(n, t, k)
    x, P, m, dw = best
    table.append(dict(design=name, v=n - t, cst_over_v=ratio, junk=s, P=P, m=m, D_over_W=dw, D_over_Wm=x, first_m=first_m))
    print(f"   {name:40s} {n-t:7d} {ratio:6.3f} {s:9.1f} {P:6d} {m:4d} {dw:6.3f} {x:9.2e} {str(first_m):>14s}")

print()
print("== (c) two-sided slack on T74 (both banks start on a line: X: U -> V^perp, Y: V -> U^perp) ==")
w, cst, D1 = 99.5556, 64.0, 35.5556
E_call = D1 / 1.6796718e-3          # E per call implied by the record (first order)
for label, Dn in (("hub's ideal lever (centre charge unchanged)", 2 * w - cst),
                  ("one extra v-drop per centre per call (my reading: T->S restore)", 2 * w - 2 * cst),
                  ("two extra v-drops per centre per call (T,S,T,S dirty pattern)", 2 * w - 3 * cst)):
    # E changes little: data climb 1 less per family (rank-1 at ln m), centres +cst per extra drop at rank 8
    extra = (Dn - (2 * w - cst)) / -cst  # number of extra drops
    dE = -w * math.log(73) + extra * cst * math.log(73 / 8)
    a = Dn / (E_call + dE)
    print(f"   {label:66s} D/call={Dn:7.2f}  a~{a:.3e}  x{a/1.6796718e-3:.2f}")

print()
print("== (d) what the 1e-1 bar needs: free coordinates per array u_eff with a = -ln(1-u/m)/ln(m/r) ==")
print(f"   {'m':>4s} {'r_eff=1':>8s} {'r_eff=4':>8s} {'r_eff=m/sqrt2':>14s}   (u_eff needed; T74 has 0.357 at m=73)")
for m in (5, 6, 8, 12, 16, 25, 41, 73):
    out = []
    for r in (1, 4, m / math.sqrt(2)):
        r = min(r, m / math.sqrt(2))
        # a = -ln(1 - u/m)/ln(m/r) = 0.1  =>  u = m (1 - (r/m)^0.1)
        out.append(m * (1 - (r / m) ** 0.1))
    print(f"   {m:4d} {out[0]:8.2f} {out[1]:8.2f} {out[2]:14.2f}")

json.dump(dict(template=table), open(os.path.join(HERE, "cap.json"), "w"), indent=1)
