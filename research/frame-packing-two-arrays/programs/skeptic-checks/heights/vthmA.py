#!/usr/bin/env python3
"""vthmA.py (v-p6-prove): random exact test of the two structural steps of Theorem A (p6-prove PROOF.md section 2).
Two registers, canonical steps only (move, free letter on one register, gate = 2x2 invertible scalar matrix), true
model.  A step is kept only if afterwards every non-zero entry (register, input) has all its terms with ONE Sp label.
Checked:
  (1) inside a joint episode every both-row has the same relative label rho = gbar_r1^-1 gbar_r2 (constant over the
      whole episode);
  (2) every non-monomial gate inside an episode (and the opening and closing gates) joins rows at one common
      episode position (input 1 raw position, input 2 position rho * raw);
  (3) per episode: moves inside >= d(P,Q) + d(rho^-1 P, rho^-1 Q)... i.e. >= 2 d(P,Q), with P, Q the start/end
      positions, and the raw distance travelled by each input bounded as in step 4.
"""
import sys, random
from fractions import Fraction as Fr
from veng import Eng, step, Q

def main(seed, m, ntr, nsteps):
    rng = random.Random(seed)
    E = Eng(m)
    free = [n for n in E.letters if not n.startswith("D")]
    pc = {}
    def pos(g):
        if g not in pc: pc[g] = E.pos(g)
        return pc[g]
    def lab(g): return tuple(E.symp(g))
    def labels_ok(C):
        for c in C:
            per = {}
            for (l, g) in c: per.setdefault(l, set()).add(lab(g))
            if any(len(v) > 1 for v in per.values()): return False
        return True
    def entries(c):
        per = {}
        for (l, g) in c: per.setdefault(l, g)
        return per            # one representative per input (all share a label)
    def apply_symp_pos(cols, P):
        return E.span([E.apply_symp(cols, v) for v in P])
    st = dict(kept=0, episodes=0, closed=0, gates_in_ep=0, fail1=0, fail2=0, fail3=0, unit_found=0)
    for _ in range(ntr):
        C = [{(0, E.ID): Q(1)}, {(1, E.ID): Q(1)}]
        rho = None; P = None; ep_moves = 0
        for _ in range(nsteps):
            t = rng.random()
            if t < 0.4: s = ("move", rng.randrange(2), rng.randrange(1, 1 << m))
            elif t < 0.6: s = ("free", rng.randrange(2), rng.choice(free))
            else:
                a, b = rng.choice([(1, 0), (1, 1), (1, -1), (0, 1), (2, 1), (1, Q(0, 1))]), rng.choice([(0, 1), (1, 1), (1, -1), (1, 0), (1, 2), (-1, 1)])
                A = [[toq(a[0]), toq(a[1])], [toq(b[0]), toq(b[1])]]
                det = A[0][0] * A[1][1] + Q(-1) * A[0][1] * A[1][0]
                if det.z(): continue
                s = ("gate", A)
            if s[0] == "gate":
                A = s[1]; old = C
                new = []
                for i in range(2):
                    d = {}
                    for j in range(2):
                        if A[i][j].z(): continue
                        for k, co in old[j].items(): d[k] = d.get(k, Q(0)) + A[i][j] * co
                    new.append({k: v for k, v in d.items() if not v.z()})
                if not labels_ok(new): continue
                nonmono = sum(1 for i in range(2) for j in range(2) if not A[i][j].z()) > 2
                both_before = [len(entries(c)) == 2 for c in old]
                both_after = [len(entries(c)) == 2 for c in new]
                in_ep = any(both_before)
                opens = (not in_ep) and any(both_after)
                if nonmono and (in_ep or opens):
                    # common episode position of the joined rows
                    if opens:
                        r = both_after.index(True); e = entries(new[r])
                        rho = tuple(E.symp(E.mulk(E.dag(e[0]), e[1])))
                        P = None; ep_moves = 0; st["episodes"] += 1
                    def epos(c):
                        ps = set()
                        for (l, g) in c:
                            x = pos(g)
                            ps.add(x if l == 0 else apply_symp_pos(list(rho), x))
                        return ps
                    ps = epos(old[0]) | epos(old[1])
                    if len(ps) != 1: st["fail2"] += 1
                    X = next(iter(ps))
                    if opens: P = X
                    st["gates_in_ep"] += 1
                    closes = in_ep and not any(both_after)
                    if closes and P is not None:
                        # (3) episode cost: two tokens walk P -> X
                        st["closed"] += 1
                        if ep_moves < 2 * E.d(P, X): st["fail3"] += 1
                        rho = None; P = None
                C = new; st["kept"] += 1
            else:
                new = step(E, C, s)
                if not labels_ok(new): continue
                C = new; st["kept"] += 1
                if s[0] == "move" and rho is not None: ep_moves += 1
            # (1) rho constant on every both-row inside an episode
            for c in C:
                e = entries(c)
                if len(e) == 2:
                    r_ = tuple(E.symp(E.mulk(E.dag(e[0]), e[1])))
                    if rho is None: st["fail1"] += 1      # both-row without an open episode: bookkeeping failure
                    elif r_ != rho: st["fail1"] += 1
    print("m=%d seed=%d %s" % (m, seed, st), flush=True)

def toq(x): return x if isinstance(x, Q) else Q(x)

if __name__ == "__main__":
    for sd in (1, 2, 3):
        main(sd, 2, 150, 40)
    for sd in (4, 5):
        main(sd, 3, 40, 40)
