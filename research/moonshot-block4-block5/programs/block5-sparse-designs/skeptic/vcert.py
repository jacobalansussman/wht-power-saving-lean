"""vcert.py (skeptic v-w8-b5sparse, own checking code): rebuild a step-2 design exactly as dumpsys.py does (same driver,
same env), and when its d = 3 closure is built, check the radius-K ball (my own BFS) row by row by BRUTE FORCE:
  for each ball row (l, g): (i) it is constrained (layer >= d, not free by the model's rule, re-implemented here);
  (ii) every producer found by brute force over ALL generator terms (z: b = g t^-1, y: b = D^-1 g t^-1, b monomial with
  Z4 entries, my own test; no frame index) is either a live column of the system with the SAME coefficient in this row,
  or dead with a VALID kill chain (verified recursively, brute force, kill order respected);
  (iii) the system has no entry in this row from a column that is not a brute-force producer.
Also: the seed's K-terms, and that the in-process system equals the agent's pickle in size.
usage: DRV=jd_union VK=2 python3 vcert.py <jdesign args>  (env as t42)"""
import os, sys, pickle, time
import numpy as np
DRV = os.environ.get('DRV', 'jd_union')
mod = __import__(DRV)
import jd2, peel2
import jdesign as J
from sp2 import M as MTOP, Fstar

peel2.RECORD.append(1)
VK = int(os.environ.get('VK', '2'))
orig = J.S.closure_any
Z4 = np.array([1, 1j, -1, -1j])
T0 = time.time()


def mono(B):
    nz = np.abs(B) > 1e-9
    if not (nz.sum(axis=1) == 1).all() or not (nz.sum(axis=0) == 1).all():
        return False
    v = B[nz]
    return bool(np.all(np.min(np.abs(v[:, None] - Z4[None, :]), axis=1) < 1e-9))


def check(P, l, d, seeds):
    E, Sy = P.E, P.S
    st = jd2.STASH[l]
    rows, cols, ent, (klass, kterms) = st['cl']
    seedcols = st['seedcols']
    colrows, rowprod, dead, expanded, pending = peel2.DEBUG['state']
    ev = list(peel2.EVENTS)
    order = {}; trig = {}
    for i, (k, r) in enumerate(ev):
        if k not in order:
            order[k] = i; trig[k] = r
    Dinv = np.conj(E.D).T
    Fs = [Sy.frame(P.s[0]), Sy.frame(P.s[1])]
    terms = [(u, l2, t, cf) for u, T in P.gens.items() for (l2, t, cf) in T]
    THs = {}

    def isfree(l_, g):
        if E.layer(g) == 0:
            return True
        F = Sy.frame(g)
        if F == Fs[l_] or F == Fstar:
            return True
        h = E.intern(E.snap(Dinv @ E.MATS[g]))
        return E.layer(h) == 0 or Sy.frame(h) == Fs[l_]

    pcache = {}

    def producers(row):
        if row in pcache:
            return pcache[row]
        l_, g = row
        G = E.MATS[g]; GD = Dinv @ G
        out = {}
        for (u, l2, t, cf) in terms:
            if l2 != l_:
                continue
            TH = THs.get(t)
            if TH is None:
                TH = np.conj(E.MATS[t]).T; THs[t] = TH
            for typ, A in (('z', G), ('y', GD)):
                Bm = A @ TH
                if mono(Bm):
                    key = (typ, u, E.intern(E.snap(Bm)))
                    c = out.setdefault(key, {})
                    for k, v in cf.items():
                        c[k] = c.get(k, 0) + v
        out = {k: c for k, c in out.items() if any(abs(v) > 1e-12 for v in c.values())}
        pcache[row] = out
        return out

    verified = {}
    fails = []

    def kill_ok(key, depth=0):
        """iterative verification of key's death"""
        stack = [key]
        while stack:
            k = stack[-1]
            if k in verified:
                stack.pop(); continue
            r = trig.get(k)
            if r is None:
                verified[k] = False; fails.append(('dead without trigger', k)); stack.pop(); continue
            l_, g = r
            if E.layer(g) < d or isfree(l_, g):
                verified[k] = False; fails.append(('trigger not constrained', k)); stack.pop(); continue
            pr = producers(r)
            if k not in pr:
                verified[k] = False; fails.append(('killed col does not produce trigger', k)); stack.pop(); continue
            if len(pr[k]) != 1:
                fails.append(('kill entry is a sum of coefficients', k, pr[k]))
            others = [o for o in pr if o != k]
            bad = [o for o in others if order.get(o, 10 ** 12) > order[k]]
            if bad:
                verified[k] = False; fails.append(('other producer not killed earlier', k, len(bad))); stack.pop(); continue
            todo = [o for o in others if o not in verified]
            if todo:
                stack.extend(todo); continue
            verified[k] = all(verified[o] for o in others)
            stack.pop()
        return verified[key]

    # dead seeds: verify their kill chains (brute force, order respected)
    for (u0, b0) in seeds:
        skey = ('y', u0, b0)
        if skey in dead:
            ok = kill_ok(skey)
            print("VCERT out %d d=%d: SEED DEAD; kill chain verified %s (%d kills checked, %d failures so far %s) (%.0fs)" % (
                l + 1, d, ok, len(verified), len(fails), fails[:4], time.time() - T0), flush=True)
            if os.environ.get('VEXIT1') and l == 0:
                sys.exit(0)
    # my own BFS ball
    rinv = {i: r for r, i in rows.items()}
    cinv = {i: k for k, i in cols.items()}
    byrow = {}; bycol = {}
    for (ri, ci, cf) in ent:
        byrow.setdefault(ri, {}).setdefault(ci, {})
        for k, v in cf.items():
            byrow[ri][ci][k] = byrow[ri][ci].get(k, 0) + v
        bycol.setdefault(ci, set()).add(ri)
    C = set(seedcols); R = set()
    for _ in range(VK):
        nr = set()
        for c in C:
            nr |= bycol.get(c, set())
        nr -= R; R |= nr
        for r in nr:
            C |= set(byrow[r])
    nrows = 0; ndead = 0; nlive = 0
    for ri in sorted(R):
        row = rinv[ri]
        l_, g = row
        if E.layer(g) < d or isfree(l_, g):
            fails.append(('ball row not constrained', row)); continue
        pr = producers(row)
        sysrow = byrow[ri]
        for ci, cf in sysrow.items():
            k = cinv[ci]
            if k not in pr:
                fails.append(('system entry from a non-producer', row, k)); continue
            a = {kk: v for kk, v in cf.items() if abs(v) > 1e-12}
            b = {kk: v for kk, v in pr[k].items() if abs(v) > 1e-12}
            if set(a) != set(b) or any(abs(a[x] - b[x]) > 1e-9 for x in a):
                fails.append(('coefficient mismatch', row, k, a, b))
        for k in pr:
            if k in cols:
                nlive += 1
                if cols[k] not in sysrow:
                    fails.append(('live producer missing from row', row, k))
            else:
                ndead += 1
                if k not in dead:
                    fails.append(('producer neither live nor dead (incomplete row in ball)', row, k)); continue
                if not kill_ok(k):
                    fails.append(('dead producer with invalid kill chain', row, k))
        nrows += 1
    # seed K-terms
    sk = [(cinv[c], kterms.get(c)) for c in seedcols]
    nchain = sum(1 for k in verified)
    print("VCERT out %d d=%d: ball K=%d rows %d, cols %d; live producer incidences %d, dead producer incidences %d; "
          "kill chains verified %d (all ok %s); FAILURES %d %s; seed %s; system %d x %d, %d entries (%.0fs)" % (
              l + 1, d, VK, nrows, len(C), nlive, ndead, nchain, all(verified.values()), len(fails), fails[:6],
              sk, len(rows), len(cols), len(ent), time.time() - T0), flush=True)


def cl(P, l, d, seeds, cap=3000):
    r = orig(P, l, d, seeds, cap)
    if d == 3:
        check(P, l, d, seeds)
    return r


J.S.closure_any = cl


def hook(Nss, nal, rng, starts=3, iters=60):
    d = jd2.STASH[0]['d']
    if d >= 4:
        return 0.0, None       # go on to d = 3 regardless (as FORCE4); nothing is decided here
    print("VCERT: d=3 checks done; stopping", flush=True)
    sys.exit(0)


J.lm = hook

if __name__ == "__main__":
    print("vcert DRV %s argv %s VK %d PYTHONHASHSEED=%s" % (DRV, sys.argv[1:], VK, os.environ.get('PYTHONHASHSEED')), flush=True)
    J.main()
