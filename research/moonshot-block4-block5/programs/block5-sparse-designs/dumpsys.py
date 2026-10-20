"""dumpsys.py (w8-b5sparse, resumed 15 Oct): rebuild a step-2 design exactly as jd_zcover.py / jd_union.py do and
PICKLE the layered systems it meets: sys_d4_<seed>.pkl (complete d = 4 closures of both outputs) and, if the numerical
d = 4 pass holds, sys_d3_<seed>.pkl (the d = 3 closure, truncated at CAP: only complete rows are a valid test).
Nothing is decided here; v4gb.py computes the exact V4 of the d = 4 system and the d = 3 local test on it.
usage: DRV=jd_zcover|jd_union python3 dumpsys.py <jdesign args>  (same env as t24 / t12b)."""
import os, sys, pickle
DRV = os.environ.get('DRV', 'jd_zcover')
mod = __import__(DRV)
import jd2
import jdesign as J


def hook(Nss, nal, rng, starts=3, iters=60):
    d = jd2.STASH[0]['d']
    pickle.dump(dict(nal=nal, d=d, drv=DRV, argv=sys.argv[1:],
                     sys=[(jd2.STASH[l]['cl'], jd2.STASH[l]['seedcols'], jd2.STASH[l]['stats'].get('cap')) for l in (0, 1)]),
                open('sys_d%d_%s.pkl' % (d, sys.argv[1]), 'wb'))
    print('     dumpsys: pickled sys_d%d_%s.pkl' % (d, sys.argv[1]), flush=True)
    if d >= 4:
        cost, al = jd2.lm2(Nss, nal, rng, starts, iters)
        if os.environ.get('FORCE4') and d == 4:      # the exact d=3 test contains the d=4 rows: go on regardless
            print('     dumpsys: FORCE4: d=4 LM residual %.3e ignored, continuing to d=3' % cost, flush=True)
            return 0.0, al
        return cost, al
    return 1.0, None


J.lm = hook

if __name__ == "__main__":
    print("dumpsys DRV %s config %s argv %s PYTHONHASHSEED=%s" % (DRV, jd2.CFG, sys.argv[1:], os.environ.get('PYTHONHASHSEED')), flush=True)
    J.main()
