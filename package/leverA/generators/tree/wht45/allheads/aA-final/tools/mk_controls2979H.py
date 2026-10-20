"""aA-final (wht45/allheads, 19 Oct 2026): the release controls of L2979H (design A), made exactly as x-final made the
#2979 ones (controls/ of lever2979/x-final), from the design-A sources:
  Heads2979HNeg.lean     x-final's Heads2979Neg.lean on aA-data's table K.DH2979H: 16,577 heads refuted, residual >= 5 refuted
  Rate2979HSharpCtl.lean e-r45m1's RateMono.lean + aA-data's generated Rate2979HSharp.lean (bodies unchanged below
                         their imports) in one file, as x-final's Rate2979SharpCtl.lean
  NegHeads.lean, NegGlob.lean   aA-data's negatives (colab-neg, Data/H2979H), unchanged
usage: python3 -I mk_controls2979H.py OUTDIR [RATE2979HSHARP.lean NEGHEADS.lean NEGGLOB.lean]
  (defaults: aA-data's files as compiled on C1; generators/REGENERATE.sh passes the regenerated ones)"""
import os, re, sys, hashlib
B = os.environ['GEN_ROOT'] + '/wht45/'
OUT = sys.argv[1]
if os.path.exists(OUT): sys.exit('REFUSE: ' + OUT + ' exists')
os.makedirs(OUT)
def rd(p): return open(B + p, encoding='utf-8').read()
strip = lambda s: re.sub(r'\A(import [^\n]*\n)+', '', s)
# 1. heads / residual control
h = rd('lever2979/x-final/controls/Heads2979Neg.lean')
h2 = h.replace('import Work.XHostT.Data.H2979.Calls', 'import Work.XHostT.Data.H2979H.Calls')
h2 = h2.replace("x-final (wht45/lever2979) Lean control, compile form: on #2979's host table (x-tables' data)",
                "aA-final (wht45/allheads) Lean control, compile form: on design A's host table (#2979 with every head kept,\naA-data's data)")
h2 = h2.replace('NOT 16,087 (one head more)', 'NOT 16,577 (one head more)')
h2 = h2.replace("`Work/XHostT/Data/H2979/Heads.lean`\nproves 16,086 and 4.  i-final's Heads3945Neg.lean with the #2979 data.",
                "`Work/XHostT/Data/H2979H/Heads.lean`\nproves 16,576 and 4.  x-final's Heads2979Neg.lean with the design-A data.")
h2 = re.sub(r'\bK\.DH2979\b', 'K.DH2979H', h2)
h2 = h2.replace('cntK rt 48 4 74 ≠ 16087', 'cntK rt 48 4 74 ≠ 16577')
assert '16087' not in h2 and '16,086' not in h2 and 'DH2979.' not in h2 and 'H2979.Calls' not in h2, h2
# 2. rate control
ctl = rd('lever2979/x-final/controls/Rate2979SharpCtl.lean')
rm = rd('lever/e-r45m1/out/bundle/Work/Mixed/RateMono.lean')
rs_old = rd('lever2979/x-inst/lean/Work/Mixed/Rate2979Sharp.lean')
assert ctl == '\n'.join(ctl.split('\n')[:3]) + '\n' + strip(rm) + strip(rs_old)     # x-final's recipe, checked
EXTRA = sys.argv[2:5] if len(sys.argv) >= 5 else None
rdx = lambda i, p: open(EXTRA[i], encoding='utf-8').read() if EXTRA else rd(p)
rs = rdx(0, 'allheads/aA-data/colab-src/T/Work/Mixed/Rate2979HSharp.lean')
assert rs.startswith('import Work.Mixed.Rate2979H\nimport Work.Mixed.RateMono\n')
r2 = ('import Work.Mixed.Rate2979H\nimport Work.BlockAccounting.Moment\n'
      "/-! aA-final: e-r45m1's Work/Mixed/RateMono.lean and aA-data's generated Work/Mixed/Rate2979HSharp.lean (both unchanged below their imports), in one file so that it compiles against the final tree without a new module. -/\n"
      + strip(rm) + strip(rs))
# 3. negatives, unchanged
ng = rdx(1, 'allheads/aA-data/colab-neg/Work/XHostT/Data/H2979H/NegHeads.lean')
gg = rdx(2, 'allheads/aA-data/colab-neg/Work/XHostT/Data/H2979H/NegGlob.lean')
for name, s in (('Heads2979HNeg.lean', h2), ('Rate2979HSharpCtl.lean', r2), ('NegHeads.lean', ng), ('NegGlob.lean', gg)):
    open(os.path.join(OUT, name), 'w', encoding='utf-8').write(s)
    print(name, hashlib.sha256(s.encode()).hexdigest()[:16], len(s.split('\n')), 'lines')
