# o-model probe P3b: drive #3815's build_fill.build() with the head-bonus hall_lp (patched copy beside build_fill.py), but
# after the LP has converged (no new Hall rows, no new columns) solve the same restricted model as a MIP (y integral),
# re-separate Hall rows on the MIP point and repeat; report LP and MIP head counts. Our code; hub files untouched.
# Args: CONSTRUCTION_DIR ORDER OUT_DIR [mip_seconds]
import sys, time, json
sys.path.insert(0, sys.argv[1])
import numpy as np, highspy
import hall_lp as H
import build_fill as BF
MIPSEC = float(sys.argv[4]) if len(sys.argv) > 4 else 900.0
def heads(lp):
    y = lp.x[lp.first_y:]; h = 0.0; frac = set()
    for (k, j), v in zip(lp.cols, y):
        if abs(v - round(v)) > 1e-6: frac.add(k)
        if k % 9 == 4 and int(H.PATS[lp.cap[k]][1][j]) == 48: h += v
    return h, frac
orig_run = H.HallLP.run
def run(self, max_rounds=200):
    r = orig_run(self, max_rounds)
    h, frac = heads(self)
    print('LP with head bonus: root %.16f  heads(LP) %.3f  fractional classes %d (span index counts %s)' % (
        r, h, len(frac), dict(sorted(__import__('collections').Counter(k % 9 for k in frac).items()))), flush=True)
    if not frac: return r
    ycols = np.arange(self.first_y, self.first_y + len(self.cols), dtype=np.int32)
    for it in range(8):
        self.h.changeColsIntegrality(len(ycols), ycols, np.array([highspy.HighsVarType.kInteger] * len(ycols)))
        self.h.setOptionValue('time_limit', MIPSEC); self.h.setOptionValue('mip_rel_gap', 1e-9)
        t0 = time.time(); self.h.run(); st = self.h.getModelStatus()
        print('MIP round %d: %s in %.0f s, obj %.6f' % (it, self.h.modelStatusToString(st), time.time() - t0, self.h.getInfo().objective_function_value), flush=True)
        sol = self.h.getSolution(); self.x = np.array(sol.col_value)
        h, frac = heads(self)
        print('   heads(MIP) %.3f  fractional classes %d  root %.16f' % (h, len(frac), H.root(self.profile())), flush=True)
        added = self.separate()
        print('   new Hall rows on the MIP point:', added, flush=True)
        if added == 0: break
        self.h.changeColsIntegrality(len(ycols), ycols, np.array([highspy.HighsVarType.kContinuous] * len(ycols)))
    return H.root(self.profile())
H.HallLP.run = run
out, r = BF.build(sys.argv[2], sys.argv[3], distinct=40)
import subprocess
res = subprocess.run([sys.executable, '-B', str(out / 'owner_audit.py')], capture_output=True, text=True)
print('owner_audit :', res.returncode, res.stdout[:600], res.stderr[-2000:], flush=True)
