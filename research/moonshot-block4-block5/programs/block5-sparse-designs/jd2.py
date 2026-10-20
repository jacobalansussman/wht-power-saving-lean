"""jd2.py (w8-b5sparse): w7-b5's jdesign.py designs (with v-w7-b5's cycle-orbit patch), decided with
  * peel2 (complete forced-zero closure; at the cap, the valid truncated system of complete rows), and
  * modsolve (exact rank test mod p = 2^61-1 at a random alpha: the generic-alpha verdict), plus
  * endsolve's dense LM over alpha (the special-alpha check) when the system has <= LMMAX columns.
Environment: CAP (columns per output, default 4000), TRUNC (1: test the truncated system at the cap), FIFO (1/0),
LMMAX (default 2500), NORM ('all' | 'seed'), OLDPEEL (1: use the original peel.py with an audit line, for comparison).
usage: as jdesign.py: jd2.py seed profA profB lam1 lam2 mode rounds nper dlist [8]"""
import os, sys, time, random
import numpy as np
import jdesign as J
import split as S
import endsolve as es
import modsolve as MS
from peel2 import peel_closure2

STASH = {}
CFG = dict(CAP=int(os.environ.get('CAP', '4000')), TRUNC=int(os.environ.get('TRUNC', '1')),
           FIFO=int(os.environ.get('FIFO', '1')), LMMAX=int(os.environ.get('LMMAX', '2500')),
           NORM=os.environ.get('NORM', 'all'), OLDPEEL=int(os.environ.get('OLDPEEL', '0')),
           PROG=int(os.environ.get('PROG', '0')), MODPMAX=float(os.environ.get('MODPMAX', '3e8')))
T0 = [time.time()]


def closure_any2(P, target_l, dmin, seeds, cap=3000):
    cap = CFG['CAP']
    if CFG['OLDPEEL']:
        import peel_old
        rows, cols, ent, kl, st = peel_old.peel_closure(P, target_l, dmin, seeds, cap=cap)
    else:
        rows, cols, ent, kl, st = peel_closure2(P, target_l, dmin, seeds, cap=cap, fifo=bool(CFG['FIFO']), verbose_every=CFG['PROG'])
    S.LAST_TRIG = [w for (w, k) in st['triggers']]
    seedcols = set(cols[('y', u0, b0)] for (u0, b0) in seeds if ('y', u0, b0) in cols)
    print("     peel%s out %d d=%d: created %d, killed %d, live %d, rows %d (all rows of live cols %s), incomplete rows %d, "
          "class live %d of %d, seeds live %d, unconstrained cols %s (class %s), requeued %s%s (%.0fs)" % (
              '_old' if CFG['OLDPEEL'] else '2', target_l + 1, dmin, st['created'], st['killed'], st['live'], st['rows'],
              st.get('rows_all'), st.get('incomplete_rows', -1), st['class_live'], st['class_total'], st['seeds_live'],
              st.get('cols_unconstrained'), st.get('class_unconstrained'), st.get('requeued'),
              ' CAP' if st.get('cap') else '', time.time() - T0[0]), flush=True)
    STASH[target_l] = dict(cl=(rows, cols, ent, kl), seedcols=seedcols, stats=st, d=dmin)
    if st.get('cap') and not CFG['TRUNC']:
        return None
    return rows, cols, ent, kl


def lm2(Nss, nal, rng, starts=3, iters=60):
    verdicts = []
    for l in (0, 1):
        s = STASH.get(l)
        if s is None:
            continue
        cl = s['cl']
        if (len(cl[0]) + 1) * (len(cl[1]) + 1) > CFG['MODPMAX']:
            print("     modp out %d d=%d: skipped (%d x %d > MODPMAX)" % (l + 1, s['d'], len(cl[0]), len(cl[1])), flush=True)
            verdicts.append(True)
            continue
        for norm in (('all', 'seed') if CFG['NORM'] == 'both' else (CFG['NORM'],)):
            t = time.time()
            ok, ra, rt, nr, nc = MS.consistent(cl, nal, seed=11 + l, norm=norm, seedcols=s['seedcols'])
            capflag = ' TRUNCATED(cap)' if s['stats'].get('cap') else ' complete'
            print("     modp out %d d=%d norm=%s%s: %d x %d, rank A %d, rank [A|e] %d -> %s at random alpha (%.0fs)" % (
                l + 1, s['d'], norm, capflag, nr, nc, ra, rt, 'CONSISTENT' if ok else 'INCONSISTENT', time.time() - t), flush=True)
            verdicts.append(ok)
    ncmax = max(Ns[0].shape[1] for Ns in Nss)
    if ncmax <= CFG['LMMAX']:
        cost, al = es.lm(Nss, nal, rng, starts=starts, iters=iters)
        print("     LM (dense, special alpha): residual %.3e" % cost, flush=True)
        if all(verdicts) and cost > 1e-12:
            print("     NOTE: modp consistent but LM residual > 0 (LM not converged?)", flush=True)
        return cost, al
    return (0.0 if all(verdicts) else 1.0), None


J.S.closure_any = closure_any2
J.lm = lm2

if __name__ == "__main__":
    print("jd2 config %s argv %s PYTHONHASHSEED=%s" % (CFG, sys.argv[1:], os.environ.get('PYTHONHASHSEED')), flush=True)
    J.main()
