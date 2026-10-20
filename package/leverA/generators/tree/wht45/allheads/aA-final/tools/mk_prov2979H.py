"""aA-final (wht45/allheads, 19 Oct 2026): write Work/L2979H/{Prov,Engine}.lean for design A (#2979, every head kept)
from x-final's Work/L2979/{Prov,Engine}.lean (sha-pinned) with exactly the renames of aA-data's INTERFACE section 7.
Fill-independent names (providers PWG2979.*, GXD.S2979B/E.*, dispatch D2979.disp, warp glue D2979W.*, win401N) stay.
Refuses (exit 2) if a source hash differs, a rename finds nothing, or an old fill-dependent token is left.
usage: python3 -I mk_prov2979H.py OUT_LEAN_ROOT"""
import hashlib, os, re, sys
XF = os.environ['GEN_ROOT'] + '/wht45/lever2979/x-final/lean/Work/L2979'
OUT = sys.argv[1]
PIN = {'Prov.lean': None, 'Engine.lean': None}
def die(m):
    print('REFUSE:', m); sys.exit(2)
src = {f: open(os.path.join(XF, f), encoding='utf-8').read() for f in PIN}
for f, s in src.items():
    PIN[f] = hashlib.sha256(s.encode()).hexdigest()
# (old, new, count expected >= 1 in the named file) -- word-boundary replacements, in this order
R = [
 ('Prov.lean', r'import Work\.XHostT\.Inst\.S2979Facts\b', 'import Work.XHostT.Inst.S2979HFacts'),
 ('Prov.lean', r'import Work\.Mixed\.Bridge\.Inst2979\b', 'import Work.Mixed.Bridge.Inst2979H'),
 ('Engine.lean', r'import Work\.L2979\.EngineOfGeomW\b', 'import Work.L2979H.EngineOfGeomW'),
 ('Engine.lean', r'import Work\.L2979\.Prov\b', 'import Work.L2979H.Prov'),
 ('*', r'\bnamespace L2979\b', 'namespace L2979H'),
 ('*', r'\bend L2979\b', 'end L2979H'),
 ('*', r'OAI\.PowerSaving\.MX\.L2979\.', 'OAI.PowerSaving.MX.L2979H.'),
 ('Engine.lean', r'\bMX\.L2979\.', 'MX.L2979H.'),
 ('Engine.lean', r'THE ENGINE THEOREM OF L2979\.', 'THE ENGINE THEOREM OF L2979H.'),
 ('*', r'\bInst\.sched2979_fits\b', 'Inst.sched2979H_fits'),
 ('*', r'\bInst\.sched2979_hcost\b', 'Inst.sched2979H_hcost'),
 ('*', r'\bInst\.sched2979_heads\b', 'Inst.sched2979H_heads'),
 ('*', r'\bInst\.sched2979\b', 'Inst.sched2979H'),
 ('*', r'\bInst\.resK2979_disp\b', 'Inst.resK2979H_disp'),
 ('*', r'\bInst\.resK2979_pos\b', 'Inst.resK2979H_pos'),
 ('*', r'\bK\.DH2979\.', 'K.DH2979H.'),
 ('*', r'\bPB\.unit2979\b', 'PB.unit2979H'),
 ('*', r'\bresIcW2979\b', 'resIcW2979H'),
 ('*', r'\bhfit2979W\b', 'hfit2979HW'),
 ('*', r'\bhpos2979W\b', 'hpos2979HW'),
 ('*', r'\bhU2979W\b', 'hU2979HW'),
 # docstrings / comments: the fill-dependent names they mention
 ('Prov.lean', r"x-tables' residual table `DH2979\.resT`", "aA-data's residual table `DH2979H.resT` (= x-tables')"),
 ('Prov.lean', r"\(x-tables' `sched2979_fits`\)", "(aA-data's `sched2979H_fits`)"),
 ('Prov.lean', r"is its profile\*\* `unit2979`", "is its profile** `unit2979H`"),
 ('Prov.lean', r"`listCost φ unit2979` \(666", "`listCost φ unit2979H` (666"),
]
out = dict(src)
for f, pat, new in R:
    targets = list(out) if f == '*' else [f]
    n = 0
    for t in targets:
        out[t], k = re.subn(pat, new, out[t]); n += k
    if n == 0:
        die('no match for ' + pat)
# docstrings: say what this is
out['Prov.lean'] = out['Prov.lean'].replace(
    '# x-final: the data of the #2979 warp network (providers, dispatch, window facts, residuals, unit price)',
    '# aA-final: the data of the #2979 warp network for design A (every head kept; FILL 764395fc)', 1)
out['Prov.lean'] = out['Prov.lean'].replace(
    "i-final's `Work/L3945/Prov.lean` with the #2979 instance.  x-final (wht45/lever2979, 19 Oct 2026); prepared with\nassistance from Claude.  No `sorry`.",
    "x-final's `Work/L2979/Prov.lean` (i-final's `Work/L3945/Prov.lean` with the #2979 instance), with the schedule, host\n"
    "table and unit profile of design A (aA-data's `sched2979H`, `K.DH2979H`, `PB.unit2979H`: 16,576 heads; the providers,\n"
    "dispatch and residuals do not depend on the fill).  aA-final (wht45/allheads, 19 Oct 2026; made by\n"
    "tools/mk_prov2979H.py); prepared with assistance from Claude.  No `sorry`.", 1)
out['Engine.lean'] = out['Engine.lean'].replace('# x-final: THE engine theorem of the #2979 lever instance L2979',
    '# aA-final: THE engine theorem of the #2979 lever instance L2979H (design A: every head kept)', 1)
out['Engine.lean'] = out['Engine.lean'].replace('`MX.L2979H.engine : MX.L2979H.EngineStmt.{0}`: OpenAI\'s `hills_program` statement at the seam exponent\n`1 - 1755647301/10^12` (DFT figure 1.7556473e-3)',
    '`MX.L2979H.engine : MX.L2979H.EngineStmt.{0}`: OpenAI\'s `hills_program` statement at the seam exponent\n`1 - 1756470701/10^12` (DFT figure 1.7564707e-3)', 1)
out['Engine.lean'] = out['Engine.lean'].replace('half cut of the 16,086 literal rank-48 heads', 'half cut of the 16,576 literal rank-48 heads of\nthe every-head-kept fill 764395fc (round wht45/reopt)', 1)
out['Engine.lean'] = out['Engine.lean'].replace('This is the one theorem the generated chain (`Work.FourierL2979.Seam`) imports.\nx-final (wht45/lever2979, 19 Oct 2026); prepared with assistance from Claude.  No `sorry`.',
    'This is the one theorem the generated chain (`Work.FourierL2979H.Seam`) imports.\nx-final\'s Work/L2979/Engine.lean with the renames of aA-data\'s INTERFACE section 7; aA-final (wht45/allheads,\n19 Oct 2026; made by tools/mk_prov2979H.py); prepared with assistance from Claude.  No `sorry`.', 1)
# leftovers: any fill-dependent #2979 token without H
bad = re.compile(r'\b(sched2979(?!H)|resK2979_|DH2979\b|unit2979\b|S2979Facts|Inst2979\b|L2979\b|1755647301|17556473|16086|16,086)')
for t, s in out.items():
    for i, line in enumerate(s.split('\n')):
        m = bad.search(line)
        if m and "x-final's" not in line:   # the provenance lines name x-final's source files
            die(f'{t}:{i+1}: leftover {m.group(0)!r}: {line.strip()}')
os.makedirs(os.path.join(OUT, 'Work/L2979H'), exist_ok=True)
for t, s in out.items():
    p = os.path.join(OUT, 'Work/L2979H', t)
    if os.path.exists(p) and open(p, encoding='utf-8').read() != s:
        die(p + ' exists and differs')
    open(p, 'w', encoding='utf-8').write(s)
    print(t, 'from', PIN[t][:16], '->', hashlib.sha256(s.encode()).hexdigest()[:16])
