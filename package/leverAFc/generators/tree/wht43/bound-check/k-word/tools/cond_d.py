"""k-word: x-mutable's condition (d) on a hub canonical word with exit frames (exact rational arithmetic).

Mutable sources (x-mutable): the sources are written in place and restored by the DERIVED end gate, which runs at the
FULL frame and sets the undone rows z = (x roles, slots) to the rows of the inverse.  With exit frames (x-exit-lean,
XExit/Hub.lean) the end gate is split: the hub's TAIL gates (one cleanup gate per exit helper, at its exit frame) run
first, and the REST is R := (T Mt)^-1 = Mti (1 + K) on the undone rows.  A gate at the FULL frame can only read a
register that is at the FULL frame, and an exit helper never climbs back up, so:

  condition (d): every SOURCE row of R reads only registers that end at the FULL frame, i.e. R[x_i][q] = 0 for every
                 exit helper q (final frame != FULL).

(In the classic engine the source rows of R are unit rows, so (d) is automatic there; with mutable sources it is new.)
Reported, besides (d):
  * the same test on the plain inverse Mti (no tail), which is what a lowering that DROPS the tail would need;
  * hRrow / hRcol of x-exit-lean for the slot rows (exit rows of R unit, exit columns of R unit);
  * the cross-check  rest-cleanup matrix == R  on the undone rows (the hub's own factorisation of R), and whether the
    rest cleanup touches only FULL-ended registers (the hub audit's 'cleanup on non-FULL end' rule).
Representation: the hub canonical word itself (gates 'add' t s num den frame: t += num/den s; 'scale' t a b frame).
The gcert lowering changes only (i) scale units (a diagonal similarity: zero patterns unchanged), (ii) target-only
gates (precondition, compensation, cut), which do not enter the undone rows, (iii) the source copies (absent here).
usage: cond_d.py CANONICAL-WORD.json.gz OUT.json"""
import sys, json, gzip
from fractions import Fraction as Q
from collections import Counter

def read(p):
    raw = open(p, 'rb').read()
    return json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)

def apply(e, rows, inverse=False):
    if e[0] == 'scale':
        _, t, a, b, *_ = e
        q = Q(a, b)
        if inverse: q = 1 / q
        rows[t] = {j: q * x for j, x in rows[t].items()}
    else:
        _, t, s, a, b, *_ = e
        assert t != s
        q = Q(a, b) * (-1 if inverse else 1)
        d = rows[t]
        for j, x in list(rows[s].items()):
            w = d.get(j, 0) + q * x
            if w: d[j] = w
            else: d.pop(j, None)

def prog_rows(n, gates):
    rows = [{r: Q(1)} for r in range(n)]
    for e in gates: apply(e, rows)
    return rows

def inv_rows(n, gates):
    """rows of the inverse of the program `gates`"""
    rows = [{r: Q(1)} for r in range(n)]
    for e in reversed(gates): apply(e, rows, inverse=True)
    return rows

def main(cp, outp):
    c = read(cp)
    v, R, h = c['v'], c['R'], c['h']; n = 2 * v + R
    F = c['frames']; FULL = 1; assert len(F[FULL]) == h
    final = c['final']
    cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'
    regs = lambda e: [e[1]] + ([e[2]] if e[0] == 'add' else [])
    main_ = c['A'] + c['cut'] + c['B']
    tail = c.get('tail', []); rest = c['cleanup']
    Y = set(range(v, 2 * v)); Z = [r for r in range(n) if r not in Y]
    exitq = sorted(r for r in range(2 * v, n) if final[r] != FULL)
    nonfull_nonexit = sorted(r for r in Z if final[r] != FULL and r not in exitq)
    # rows of Mti and of R = (T Mt)^-1  (program: tail^-1 then main^-1)
    Mti = inv_rows(n, main_)
    Rr = inv_rows(n, main_ + tail)
    assert all(not (set(Rr[r]) & Y) for r in Z), 'an undone row of R reads a target'
    assert all(not (set(Mti[r]) & Y) for r in Z), 'an undone row of Mti reads a target'
    ex = set(exitq)
    def viol(rows, rowset):
        out = []
        for r in rowset:
            bad = sorted(j for j in rows[r] if j in ex and j != r)
            if bad: out.append([r, bad])
        return out
    xs = list(range(v))
    d_R = viol(Rr, xs)
    d_Mti = viol(Mti, xs)
    fullslots = [q for q in range(2 * v, n) if q not in ex]
    hRcol_slots = viol(Rr, fullslots)
    hRrow = [q for q in exitq if Rr[q] != {q: Q(1)}]
    hRcol = [q for q in exitq if any(Rr[r].get(q) for r in Z if r != q)]
    # the hub's rest cleanup as a matrix on the undone rows, and its registers' final frames
    Rest = prog_rows(n, rest)
    rest_eq = all(Rest[r] == Rr[r] for r in Z)
    rest_regs_nonfull = sorted({r for e in rest for r in regs(e) if final[r] != FULL})
    # partners of the source rows (who the derived end gate reads for x)
    part = Counter()
    for i in xs:
        for j in Rr[i]:
            if j != i: part[cls(j) + (':exit' if j in ex else ':full')] += 1
    tail_cols = Counter(cls(e[2]) + (':exit' if e[2] in ex else ':full') for e in tail if e[0] == 'add')
    tail_rows_exit = all(e[1] in ex for e in tail)
    xrows_unit_Mti = sum(1 for i in xs if Mti[i] == {i: Q(1)})
    res = dict(word=cp, v=v, R=R, exit_helpers=len(exitq),
               exit_dims=dict(sorted(Counter(len(F[final[q]]) for q in exitq).items())),
               nonFULL_nonexit_undone_registers=nonfull_nonexit,
               tail_gates=len(tail), tail_rows_are_exit_helpers=tail_rows_exit, tail_partner_kinds=dict(tail_cols),
               source_rows_unit_in_Mti=xrows_unit_Mti, source_rows_mutated=v - xrows_unit_Mti,
               condition_d_holds=not d_R, condition_d_violations=d_R,
               condition_d_without_tail_holds=not d_Mti, condition_d_without_tail_violations=d_Mti[:40],
               condition_d_without_tail_violation_count=len(d_Mti),
               source_row_partner_kinds=dict(part),
               hRrow_exit_rows_unit=not hRrow, hRrow_failures=hRrow[:20],
               hRcol_exit_columns_unit=not hRcol, hRcol_failures=hRcol[:20],
               fullslot_rows_reading_exit=len(hRcol_slots),
               rest_cleanup_equals_R_on_undone_rows=rest_eq,
               rest_cleanup_touches_nonFULL_registers=rest_regs_nonfull[:20])
    print(json.dumps(res, indent=1))
    json.dump(res, open(outp, 'w'), indent=1)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
