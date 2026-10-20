"""afc-final (wht45/afc, 20 Oct 2026): the release controls of L2979Fc (A-fc: #2979's words, the re-solved fill
u1-2979sep, the full cut), made exactly as aA-final made design A's (allheads/aA-final/tools/mk_controls2979H.py),
from the A-fc sources:
  Heads2979FcNeg.lean     aA-final's Heads2979HNeg.lean on afc-data's table K.DH2979Fc: 16,577 heads refuted,
                          residual >= 5 refuted (A-fc's residuals are x-tables', smallest width 4)
  Rate2979FcSharpCtl.lean e-r45m1's RateMono.lean + afc-cut's generated Rate2979FcSharp.lean (bodies unchanged below
                          their imports) in one file, as aA-final's Rate2979HSharpCtl.lean
  Ctl2979Fc.lean          afc-cut's full-cut profile-identity controls (lean-ctl/Work/Mixed/Bridge/Ctl2979Fc.lean), unchanged
  NegHeads.lean, NegGlob.lean, NegF1.lean   afc-data's negatives (colab-neg, Data/H2979Fc), unchanged
  perturbed/ (NegSome119, NegNone119, NegSome46, NegNone46, NegSome334, NegNone334)   afc-data's, unchanged
usage: python3 -I mk_controls2979Fc.py OUTDIR [RATE2979FCSHARP.lean NEGDIR(colab-neg root)]
  (defaults: afc-cut's and afc-data's files as compiled on C4; a regeneration passes the regenerated ones)"""
import os, re, sys, hashlib
B = os.environ['GEN_ROOT'] + '/wht45/'
OUT = sys.argv[1]
if os.path.exists(OUT): sys.exit('REFUSE: ' + OUT + ' exists')
def rd(p): return open(B + p, encoding='utf-8').read()
def sha(s): return hashlib.sha256(s.encode()).hexdigest()
strip = lambda s: re.sub(r'\A(import [^\n]*\n)+', '', s)
PIN = {  # sources, sha-pinned
    'allheads/aA-final/controls/Heads2979HNeg.lean': '3b96dcd49a3b9459',
    'allheads/aA-final/controls/Rate2979HSharpCtl.lean': '6c18ff66417445b5',
    'lever/e-r45m1/out/bundle/Work/Mixed/RateMono.lean': 'df1f0ddc4d8b8201',
}
for p, h in PIN.items():
    if not sha(rd(p)).startswith(h): sys.exit('REFUSE: pinned source changed: ' + p)
# 1. heads / residual control (aA-final's, renamed to A-fc's table)
h = rd('allheads/aA-final/controls/Heads2979HNeg.lean')
h2 = h.replace('import Work.XHostT.Data.H2979H.Calls', 'import Work.XHostT.Data.H2979Fc.Calls')
h2 = h2.replace("aA-final (wht45/allheads) Lean control, compile form: on design A's host table (#2979 with every head kept,\naA-data's data)",
                "afc-final (wht45/afc) Lean control, compile form: on A-fc's host table (#2979's words with the re-solved fill\nu1-2979sep, afc-data's data)")
h2 = h2.replace("`Work/XHostT/Data/H2979H/Heads.lean`\nproves 16,576 and 4.  x-final's Heads2979Neg.lean with the design-A data.",
                "`Work/XHostT/Data/H2979Fc/Heads.lean`\nproves 16,576 and 4.  aA-final's Heads2979HNeg.lean with the A-fc data.")
h2 = re.sub(r'\bK\.DH2979H\b', 'K.DH2979Fc', h2)
assert h2 != h and 'H2979H' not in h2 and 'design A' not in h2 and 'aA-final (wht45' not in h2, h2
assert 'cntK rt 48 4 74 ≠ 16577' in h2 and 'import Work.XHostT.Data.H2979Fc.Calls' in h2
# 2. rate control (aA-final's recipe, checked on its own output first)
ctlA = rd('allheads/aA-final/controls/Rate2979HSharpCtl.lean')
rm = rd('lever/e-r45m1/out/bundle/Work/Mixed/RateMono.lean')
rsA = rd('allheads/aA-data/colab-src/T/Work/Mixed/Rate2979HSharp.lean')
hdrA = ('import Work.Mixed.Rate2979H\nimport Work.BlockAccounting.Moment\n'
        "/-! aA-final: e-r45m1's Work/Mixed/RateMono.lean and aA-data's generated Work/Mixed/Rate2979HSharp.lean (both unchanged below their imports), in one file so that it compiles against the final tree without a new module. -/\n")
assert ctlA == hdrA + strip(rm) + strip(rsA)      # aA-final's recipe, checked
EXTRA = sys.argv[2:4] if len(sys.argv) >= 4 else None
rs = open(EXTRA[0], encoding='utf-8').read() if EXTRA else rd('afc/afc-cut/colab-src/eng/Work/Mixed/Rate2979FcSharp.lean')
assert rs.startswith('import Work.Mixed.Rate2979Fc\nimport Work.Mixed.RateMono\n')
assert 'rate2979Fc_fails' in rs and '17887046' in rs
r2 = ('import Work.Mixed.Rate2979Fc\nimport Work.BlockAccounting.Moment\n'
      "/-! afc-final: e-r45m1's Work/Mixed/RateMono.lean and afc-cut's generated Work/Mixed/Rate2979FcSharp.lean (both unchanged below their imports), in one file so that it compiles against the final tree without a new module. -/\n"
      + strip(rm) + strip(rs))
# 3. afc-cut's profile-identity controls, unchanged
c3 = rd('afc/afc-cut/lean-ctl/Work/Mixed/Bridge/Ctl2979Fc.lean')
assert c3 == rd('afc/afc-cut/colab-src/ctl/Work/Mixed/Bridge/Ctl2979Fc.lean')
# 4. afc-data's negatives, unchanged
ND = EXTRA[1] if EXTRA else B + 'afc/afc-data/colab-neg'
rn = lambda p: open(os.path.join(ND, p), encoding='utf-8').read()
files = [('Heads2979FcNeg.lean', h2), ('Rate2979FcSharpCtl.lean', r2), ('Ctl2979Fc.lean', c3),
         ('NegHeads.lean', rn('Work/XHostT/Data/H2979Fc/NegHeads.lean')),
         ('NegGlob.lean', rn('Work/XHostT/Data/H2979Fc/NegGlob.lean')),
         ('NegF1.lean', rn('Work/XHostT/Data/H2979Fc/NegF1.lean'))]
for k in ('119',):
    for kind in ('Some', 'None'):
        files.append(('perturbed/Neg%s%s.lean' % (kind, k), rn('Work/XHostT/Data/H2979Fcnload/Neg%s%s.lean' % (kind, k))))
for k in ('46', '334'):
    for kind in ('Some', 'None'):
        files.append(('perturbed/Neg%s%s.lean' % (kind, k), rn('Work/XHostT/Data/H2979Fcnwin/Neg%s%s.lean' % (kind, k))))
assert '16577' in files[3][1] and 'f1_cnt' in files[5][1] or 'cntK rt 48 6 74' in files[5][1], files[5][1][:400]
os.makedirs(OUT + '/perturbed')
for name, s in files:
    open(os.path.join(OUT, name), 'w', encoding='utf-8').write(s)
    print(name, sha(s)[:16], len(s.split('\n')), 'lines')
