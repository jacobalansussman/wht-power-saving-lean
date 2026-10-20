"""famcheck.py (w7-b5): for an orbit design that passes d = 4, list the live columns of the d = 4 peel closure and the
branch through which each one's K-term arises (coefficient 1 = branch a, mu = 0.7 = branch b), to test whether the live
closure is ONE class family (Corollary O applies) or two.  usage: famcheck.py seed profA profB lam1 lam2 mode dec"""
import sys, random
import jdesign as J
import split as S
import peel as PL
import endsolve as es, designer as dz

def main():
    a = sys.argv[1:]
    seed = int(a[0]); pp = lambda s: None if s == '0' else (int(s[0]), int(s[1]))
    profA, profB = pp(a[1]), pp(a[2]); lam1 = None if a[3] == '0' else int(a[3]); lam2 = None if a[4] == '0' else int(a[4])
    mode = a[5]; dec8 = len(a) > 6 and a[6] == '8'
    # rerun the design exactly as jdesign.run does, but stop at d = 4 and inspect
    orig = S.closure_any
    info = {}
    def spy(P, target_l, dmin, seeds, cap=3000):
        if dmin == 4:
            rows, cols, ent, kl, st = PL.peel_closure(P, target_l, dmin, seeds, cap=cap)
            klass, kterms = kl
            fam = {}
            for c, terms in kterms.items():
                for ((l2, g3), cf) in terms:
                    v = list(cf.values())[0] if cf else None
                    key = ('a' if abs(v - 1.0) < 1e-9 else 'b' if abs(v - 0.7) < 1e-9 else str(v))
                    fam[key] = fam.get(key, 0) + 1
            gens = {}
            for k in cols:
                gens[(k[0], k[1])] = gens.get((k[0], k[1]), 0) + 1
            info[target_l] = dict(live=len(cols), klass=len(klass), Kterm_branches=fam, live_by_kind=gens)
            print("famcheck out %d: %s" % (target_l + 1, info[target_l]), flush=True)
        return orig(P, target_l, dmin, seeds, cap)
    S.closure_any = spy
    J.S.closure_any = spy
    r = J.run(seed, profA, profB, lam1, lam2, mode, 0, 0, [5, 4], dec8)
    print("result", r)

main()
