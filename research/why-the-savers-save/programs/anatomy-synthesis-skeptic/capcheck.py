"""v-an-synth: re-derive an-synth's cap numbers and proposal-1 thresholds from its stated formulas (own code).
Run: python3 -I calc/capcheck.py"""
import math, json
out = {}
# an-synth's template model (NOTE sec. 3, cap.py (b)): m = 1+(k-1)(2P+1), D/W = 1 - (2P+1) cst/(2P v), u = 1
designs = [('P112 (E8)', 112, 64, 9), ('E7 (56 ports)', 56, 42, 7)]
for h in (12, 16, 20, 25, 30):
    designs.append((f'Alon h={h}', math.comb(h, 3), (h + 1) * h, h + 1))
rows = []
for name, v, cst, k in designs:
    for P in range(1, 9):
        m = 1 + (k - 1) * (2 * P + 1)
        dw = 1 - (2 * P + 1) * cst / (2 * P * v)
        if dw <= 0: continue
        x = dw / m
        rows.append(dict(design=name, P=P, m=m, D_over_W=dw, D_over_Wm=x,
                         a_if_lbar_ln_m_over_9=x / math.log(m / 9),
                         a_if_lbar_ln_m_over_3p94=x / math.log(m / 3.9448),
                         a_exact_rank_r9=-math.log(1 - dw / m) / math.log(m / 9)))
best9 = max(rows, key=lambda r: r['a_if_lbar_ln_m_over_9'])
best4 = max(rows, key=lambda r: r['a_if_lbar_ln_m_over_3p94'])
bestx = max(rows, key=lambda r: r['D_over_Wm'])
print('max D/(Wm) over designs x P:', bestx['design'], bestx['P'], bestx['m'], f"{bestx['D_over_Wm']:.3e}")
print('an-synth: "with lbar >= ln(m/9), that gives a <~ 3e-3".  Its own relaxation gives:')
print(f"   max a with lbar = ln(m/9):    {best9['a_if_lbar_ln_m_over_9']:.3e}  at {best9['design']} P={best9['P']} m={best9['m']}")
print(f"   max a with lbar = ln(m/3.94): {best4['a_if_lbar_ln_m_over_3p94']:.3e}  at {best4['design']} P={best4['P']} m={best4['m']}")
for r in rows:
    if r['design'].startswith('P112') and r['P'] <= 4:
        print(f"   P112 P={r['P']} m={r['m']}: D/W={r['D_over_W']:.4f} D/(Wm)={r['D_over_Wm']:.3e}  a(ln m/9)={r['a_if_lbar_ln_m_over_9']:.3e}  a(ln m/3.94)={r['a_if_lbar_ln_m_over_3p94']:.3e}")
out['cap'] = dict(best_lbar9=best9, best_lbar394=best4, best_x=bestx)

# 1e-1 bar: u = m (1 - (r/m)^0.1); an-synth's own column r = m/sqrt2
print('\n1e-1 bar, u_eff needed at r = m/sqrt2 (an-synth table (d) third column) and a at u_eff = 1:')
for m in (25, 29, 30, 41, 73):
    u = m * (1 - (1 / math.sqrt(2)) ** 0.1)
    a1 = -math.log(1 - 1 / m) / math.log(math.sqrt(2))
    print(f'   m={m}: u_needed={u:.3f}   a(u=1, r=m/sqrt2)={a1:.4f}')
out['bar_r_m_sqrt2'] = {m: m * (1 - (1 / math.sqrt(2)) ** 0.1) for m in (25, 41, 73)}

# proposal 1: break-even extra charge per call (an-synth cap.py (c) model)
w, cst, D1, arec = 99.5556, 64.0, 35.5556, 1.6796718e-3
E_call = D1 / arec
def a_two(extra):
    Dn = 2 * w - cst - extra
    dE = -w * math.log(73) + (extra / cst) * cst * math.log(73 / 8)
    return Dn / (E_call + dE)
lo, hi = 0.0, 200.0
for _ in range(80):
    mid = (lo + hi) / 2
    if a_two(mid) > arec: lo = mid
    else: hi = mid
print(f'\nproposal 1: break-even extra charge per call (a = record) = {lo:.1f}; kill mark used by an-synth = 128; pass mark = 64')
for e in (0, 64, 99.6, 110, 128):
    print(f'   extra {e:6.1f}/call -> a = {a_two(e):.3e}  x{a_two(e)/arec:.2f}')
# same, with the m=74 cost (an-build: s is a new coordinate) ignored here -> optimistic
out['p1_break_even_extra_per_call'] = lo
out['p1_a_at'] = {str(e): a_two(e) for e in (0, 64, 99.6, 110, 128)}
# proposal 3 check
a3 = -math.log(1 - 2 / 16) / math.log(16 / 4)
print(f'\nproposal 3: u=2, m=16, r=4: a = {a3:.4f} (< 0.1)')
out['p3'] = a3
json.dump(out, open('calc/capcheck.json', 'w'), indent=1, default=str)
