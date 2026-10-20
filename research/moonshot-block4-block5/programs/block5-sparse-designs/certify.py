"""certify.py (w8-b5sparse): INDEPENDENT verification of peel2's forced-zero certificates at d = 3.

peel2 records every kill as (column key, trigger row).  This script re-checks, with its own brute-force producer
enumeration (every term t of every generator u: z-producer b = g t^-1, y-producer b = D^-1 g t^-1, kept iff b is a
monomial Clifford, i.e. a letter of B), that for each kill event in order:
  (a) the trigger row is a CONSTRAINED row (layer >= d and not free),
  (b) every producer of the trigger row other than the killed column was killed EARLIER, and
  (c) the killed column's entry on the trigger row is not identically zero (sum of its generator coefficients there;
      a single alpha_k counts as nonzero, since every designed stream has alpha_k != 0).
Then the death of a seed / class column is an exact statement for every alpha with all alpha_k != 0.
usage: as jd2.py; certifies the d = 3 closures of both outputs (dlist 543)."""
import sys
import numpy as np
import jd2, peel2
import jdesign as J
from sp2 import M, Fstar

peel2.RECORD.append(1)
orig = J.S.closure_any
RES = []


def certify(P, l, d, seeds):
    E, Sy = P.E, P.S
    Dm = E.D; Dinv = np.conj(Dm).T
    Fs = [Sy.frame(P.s[0]), Sy.frame(P.s[1])]
    terms = [(u, l2, g2, cf) for u, T in P.gens.items() for (l2, g2, cf) in T]

    def free(l_, g):
        if E.layer(g) == 0:
            return True
        F = Sy.frame(g)
        if F == Fs[l_] or F == Fstar:
            return True
        h = E.intern(E.snap(Dinv @ E.MATS[g]))
        return E.layer(h) == 0 or Sy.frame(h) == Fs[l_]

    def producers(row):
        l_, g = row
        G = E.MATS[g]; out = {}
        for (u, l2, t, cf) in terms:
            if l2 != l_:
                continue
            TH = np.conj(E.MATS[t]).T
            for typ, pre in (('z', np.eye(E.n)), ('y', Dinv)):
                bm = E.snap(pre @ G @ TH)
                if E.is_mono(bm):
                    key = (typ, u, E.intern(bm))
                    c = out.setdefault(key, {})
                    for k, v in cf.items():
                        c[k] = c.get(k, 0) + v
        return {k: c for k, c in out.items() if any(abs(v) > 1e-12 for v in c.values())}

    killed = set(); bad = []; nonmono = 0
    ev = list(peel2.EVENTS)
    for i, (k, r) in enumerate(ev):
        if r is None:
            bad.append((i, 'no trigger')); killed.add(k); continue
        l_, g = r
        if E.layer(g) < d or free(l_, g):
            bad.append((i, 'trigger row not constrained')); killed.add(k); continue
        pr = producers(r)
        if k not in pr:
            bad.append((i, 'killed column does not produce its trigger row (or zero entry)')); killed.add(k); continue
        if len(pr[k]) > 1 or 0 in pr[k] and len(pr[k]) > 1:
            nonmono += 1
        others = set(pr) - {k}
        if not others <= killed:
            bad.append((i, 'trigger row has a live producer outside the chain: %d' % len(others - killed)))
        killed.add(k)
    # trace the death of the first seed: the DAG of kills it rests on
    trig = {k: r for (k, r) in ev}
    order = {k: i for i, (k, r) in enumerate(ev)}
    seedkeys0 = [('y', u0, b0) for (u0, b0) in seeds]
    for sk in seedkeys0[:1]:
        if sk not in trig:
            continue
        need = [sk]; seen = set(); lines = []
        while need and len(lines) < 40:
            k = need.pop()
            if k in seen:
                continue
            seen.add(k)
            r = trig.get(k)
            if r is None:
                continue
            l_, g = r
            pr = producers(r)
            others = sorted((set(pr) - {k}), key=lambda x: order.get(x, -1))
            # what kind of row: which generator terms of k land on it
            kind = []
            typ, u, b = k
            pre = E.D @ E.MATS[b] if typ == 'y' else E.MATS[b]
            for (l2, t, cf) in P.gens[u]:
                if l2 == l_ and E.intern(E.snap(pre @ E.MATS[t])) == g:
                    kind.append('main' if cf == {0: 1.0} else 'stream a%s' % list(cf)[0])
            lines.append("      kill #%d %s %s by row (array %d, layer %d) [via %s]; other producers %s" % (
                order[k], typ, u, l_ + 1, E.layer(g), ','.join(kind), [(o[0], o[1], '#%d' % order.get(o, -1)) for o in others]))
            need += others
        print("   CHAIN of seed %s (%d steps shown):" % (sk[:2], len(lines)))
        for x in lines:
            print(x)
    seedkeys = [('y', u0, b0) for (u0, b0) in seeds]
    sd = [s in killed for s in seedkeys]
    print("   CERTIFY d=%d out %d: %d kill events, %d failures %s, entries that are sums of >1 coefficient: %d; seeds killed %s"
          % (d, l + 1, len(ev), len(bad), bad[:5], nonmono, sd), flush=True)
    RES.append((d, l, len(ev), len(bad), sd))


def cl(P, l, d, seeds, cap=3000):
    r = orig(P, l, d, seeds, cap)
    if d == 3:
        certify(P, l, d, seeds)
    return r


J.S.closure_any = cl
if __name__ == "__main__":
    J.main()
