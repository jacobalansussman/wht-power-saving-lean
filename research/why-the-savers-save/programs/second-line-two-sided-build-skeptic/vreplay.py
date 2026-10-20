#!/usr/bin/env python3
"""v-an-build (skeptic of an-build): my own exact replay of a gcert/1 word and the two s-barriers.

Own code, stdlib only (Fractions).  Written from the gcert/1 field meanings in k-lower/toolkit/tools/gx/gx.py
(gate_regs: "out" = targets t += (a/b) * src, src = g[2], extra regs g[4] share the frame; "in" = g[2] += sum (a/b) s;
"neg" = negate g[2]) and from the Lean word in Work/GCert/Chain/InvThm.lean:
  forward : g1 Mt at start (y -= K s) | phase A | c += Cc s; copies down; y += Jr c | phase B | ... | g8 Mti at end
  backward: (g1i Mt)^T at start | ... | (h8 Mt)^T at end   (GRoute.gate: every nonzero off-diagonal entry (i, j)
            needs IDENTICAL matrices on roles i and j, Work/GFrame/Engine/Route.lean)
  g1 M = 1 - PY M PS, g1i M = 1 + PY M PS, g8 Mi = PX + PY + PS Mi, h8 M = PX + PY + PS M   (Carrier/Live.lean)
So the nonzero pairs of the start gate are {(y_T, q): Mt[y_T, q] != 0} in BOTH words, and of the end gate
{(q, k): Mti[q, k] != 0} (forward) and {(k, q): Mt[q, k] != 0} (backward).

usage: python3 -I vreplay.py WORD.json.gz OUT.json
"""
import gzip, json, sys
from fractions import Fraction as F


def load(p):
    return json.load(gzip.open(p, 'rt') if p.endswith('.gz') else open(p))


def gate_parts(g):
    """(registers in the gate, list of elementary adds (target, source, coeff), negated register or None)"""
    if g[0] == 'out':
        regs = [g[2]] + [t for t, _, _ in g[3]] + list(g[4])
        return regs, [(t, g[2], F(a, b)) for t, a, b in g[3]], None
    if g[0] == 'in':
        regs = [g[2]] + [s for s, _, _ in g[3]]
        return regs, [(g[2], s, F(a, b)) for s, a, b in g[3]], None
    if g[0] == 'neg':
        return [g[2]], [], g[2]
    raise ValueError(g[0])


class Rows:
    """row[r] = the current content of register r as a sparse linear form in the initial registers."""
    def __init__(self, n):
        self.row = [{r: F(1)} for r in range(n)]

    def add(self, t, s, c):          # r_t += c * r_s   (s != t)
        assert s != t
        rt = self.row[t]
        for k, val in self.row[s].items():
            nv = rt.get(k, 0) + c * val
            if nv:
                rt[k] = nv
            else:
                rt.pop(k, None)

    def addform(self, t, form, c):
        rt = self.row[t]
        for k, val in form.items():
            nv = rt.get(k, 0) + c * val
            if nv:
                rt[k] = nv
            else:
                rt.pop(k, None)

    def neg(self, r):
        self.row[r] = {k: -v for k, v in self.row[r].items()}


def forward(d):
    v, R = d['v'], d['R']; n = 2 * v + R
    X = Rows(n)
    for g in d['A']:
        _, adds, ng = gate_parts(g)
        for t, s, c in adds: X.add(t, s, c)
        if ng is not None: X.neg(ng)
    copies = {k: dict(X.row[src]) for k, src, _ in d['ret']}
    for T, ent in enumerate(d['scat']['table']):
        for k, a, b in ent:
            X.addform(v + T, copies[k], F(a, b))
    for g in d['B']:
        _, adds, ng = gate_parts(g)
        for t, s, c in adds: X.add(t, s, c)
        if ng is not None: X.neg(ng)
    return X.row


def inverse(d):
    """rows of M^-1: undo B in reverse, undo the scatter (copies re-read at that point), undo A in reverse."""
    v, R = d['v'], d['R']; n = 2 * v + R
    X = Rows(n)
    for g in reversed(d['B']):
        _, adds, ng = gate_parts(g)
        if ng is not None: X.neg(ng)
        for t, s, c in reversed(adds): X.add(t, s, -c)
    copies = {k: dict(X.row[src]) for k, src, _ in d['ret']}
    for T, ent in enumerate(d['scat']['table']):
        for k, a, b in ent:
            X.addform(v + T, copies[k], -F(a, b))
    for g in reversed(d['A']):
        _, adds, ng = gate_parts(g)
        if ng is not None: X.neg(ng)
        for t, s, c in reversed(adds): X.add(t, s, -c)
    return X.row


def rank(vecs):
    """exact rank of sparse rational vectors (dicts), plain Gaussian elimination."""
    piv = {}
    rk = 0
    for vec in vecs:
        vec = {k: F(x) for k, x in vec.items() if x}
        while vec:
            p = min(vec)
            if p not in piv:
                piv[p] = vec; rk += 1
                break
            pv = piv[p]; c = vec[p] / pv[p]
            for k, val in pv.items():
                nv = vec.get(k, 0) - c * val
                if nv: vec[k] = nv
                else: vec.pop(k, None)
    return rk


def main():
    d = load(sys.argv[1])
    v, R = d['v'], d['R']; n = 2 * v + R
    isx = lambda r: r < v
    isy = lambda r: v <= r < 2 * v
    iss = lambda r: r >= 2 * v
    M = forward(d)
    Mi = inverse(d)
    out = dict(word=sys.argv[1], v=v, R=R, gates_A=len(d['A']), gates_B=len(d['B']))
    # product check M * Mi = 1 on a sample of rows (exact)
    bad = 0
    for r in range(0, n, 7):
        acc = {}
        for k, a in M[r].items():
            for j, b in Mi[k].items():
                acc[j] = acc.get(j, 0) + a * b
        acc = {j: x for j, x in acc.items() if x}
        if acc != {r: F(1)}: bad += 1
    out['M_times_Minv_identity_on_sampled_rows'] = (bad == 0)
    # Lean-type facts: x untouched by the main phase; y = y + K x + (slot terms); slots restored?
    out['x_rows_identity'] = all(M[r] == {r: F(1)} for r in range(v))
    out['y_dep_on_x'] = sum(1 for T in range(v) if any(isx(k) for k in M[v + T]))
    out['y_dep_on_other_y'] = sum(1 for T in range(v) if any(isy(k) and k != v + T for k in M[v + T]))
    out['y_dep_on_slots'] = sum(1 for T in range(v) if any(iss(k) for k in M[v + T]))
    out['slots_restored_by_main_phase'] = sum(1 for q in range(2 * v, n) if M[q] == {q: F(1)})
    # ---- start barrier (same pairs forward and backward) ----
    dirt = sorted(set(q for T in range(v) for q in M[v + T] if iss(q)))
    out['dirt_slots_joined_with_y_at_start'] = len(dirt)
    dset = set(dirt)
    metx_A = set()
    for g in d['A']:
        regs, _, _ = gate_parts(g)
        if any(isx(r) for r in regs):
            metx_A.update(r for r in regs if r in dset)
    out['dirt_slots_meeting_x_in_A'] = len(metx_A)
    Kx = sorted(metx_A)
    out['rank_start_map_y_to_Kx'] = rank([{T: M[v + T][q] for T in range(v) if q in M[v + T]} for q in Kx])
    # also the full start map onto ALL dirt slots, for reference
    out['rank_start_map_y_to_all_dirt'] = rank([{T: M[v + T][q] for T in range(v) if q in M[v + T]} for q in dirt])
    # ---- s propagation in time order, from the y's only (no dirt gate, no descents) ----
    hold = [False] * n
    for r in range(v, 2 * v): hold[r] = True
    conflict_x_gates = {'A': 0, 'scatter': 0, 'B': 0}
    for ph in ('A', 'B'):
        if ph == 'B':
            # the scatter: the centre copies c_k (separate roles) read slot ret[k] and are added into y.
            # c += Cc s happens with the slot at its frame; y += Jr c joins c with y: c holds s, the slot does not.
            pass
        for g in d[ph]:
            regs, _, _ = gate_parts(g)
            if any(hold[r] for r in regs):
                if any(isx(r) for r in regs): conflict_x_gates[ph] += 1
                for r in regs:
                    if not isx(r): hold[r] = True
    out['gates_joining_x_with_an_s_holder_in_phase'] = conflict_x_gates
    forced = [q for q in range(2 * v, n) if hold[q]]
    out['slots_forced_to_hold_s_by_y_gates'] = len(forced)
    # direct x-y gates
    dxy = {ph: sum(1 for g in d[ph] if any(isx(r) for r in gate_parts(g)[0]) and any(isy(r) for r in gate_parts(g)[0]))
           for ph in ('A', 'B')}
    out['direct_x_y_gates'] = dxy
    # ---- end barrier, backward word: (h8 Mt)^T pairs (x_k, q) with Mt[q, x_k] != 0 ----
    Fb = [q for q in forced if any(isx(k) for k in M[q])]
    out['end_bwd_forced_slots_with_x_dependence'] = len(Fb)
    out['rank_end_bwd'] = rank([{k: x for k, x in M[q].items() if isx(k)} for q in Fb])
    # slots joined to BOTH x and y at the backward end gate (one gate, identical frames => x frame = y frame)
    out['end_bwd_slots_with_x_and_y_entries'] = sum(1 for q in range(2 * v, n)
                                                    if any(isx(k) for k in M[q]) and any(isy(k) for k in M[q]))
    # ---- end barrier, forward word: g8 Mti pairs (q, x_k) with Mti[q, x_k] != 0 ----
    Ff = [q for q in forced if any(isx(k) for k in Mi[q])]
    out['end_fwd_forced_slots_with_x_dependence'] = len(Ff)
    out['rank_end_fwd'] = rank([{k: x for k, x in Mi[q].items() if isx(k)} for q in Ff])
    out['end_fwd_slots_with_x_and_y_entries'] = sum(1 for q in range(2 * v, n)
                                                    if any(isx(k) for k in Mi[q]) and any(isy(k) for k in Mi[q]))
    # ---- the y -> x information a GB call must move (x -= B y): rank of the shear block ----
    # forward shear block K = M[y, x]; GB is its transposed inverse, rank equal.
    out['rank_shear_block_K'] = rank([{k: x for k, x in M[v + T].items() if isx(k)} for T in range(v)])
    for k, val in out.items():
        print('%-48s %s' % (k, val))
    json.dump(out, open(sys.argv[2], 'w'), indent=1)


if __name__ == '__main__':
    main()
