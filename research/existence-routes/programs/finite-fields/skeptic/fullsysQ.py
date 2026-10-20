# fullsysQ.py -- v-ex2-modp: nag-rerun's fullsys.py patched to EXACT rational arithmetic (no prime), extras off.
"""fullsys.py -- nag-rerun: exact decision of the FULL exact-structure system of a (sub-)skeleton, in any zero-pattern
chart, modulo the two primes.  Extends v-nag's unitelim.py (same closure, same constant-pivot idea) with:
  * --keep 'move:j,j;move:j'  drop every summand not listed (the relaxation forces those coefficients to 0);
  * --mode u : EVERY remaining coefficient is a variable (no normalisation; covers every zero pattern at once)
    --mode n : the first remaining summand of each move has coefficient 1 (a chart; valid only together with the
               charts where that summand is 0, which are other --keep specs);
  * --pins 'k=v,...'  fix ORIGINAL alpha k (0-based, recipe numbering) to the rational v;
  * --extra GBFILE    append the alpha-only elements of the un-normalised RELAXATION's Groebner basis (msolve output,
               same prime, original numbering a1..aN): they lie in the relaxation ideal, which every unit satisfies
               (Lemma 1), so adding them is sound; dropped / pinned alphas are substituted;
  * elimination and export are done mod p directly (fast), one prime at a time, using ONLY pivots whose coefficient
    is a nonzero constant (an equivalence for every alpha); the remaining core goes to msolve -g 2.
Usage: fullsys.py recipe.json OUTPREFIX [--keep S] [--mode u|n] [--pins S] [--extra1 F] [--extra2 F] [--arrays 0,1]
                  [--cap N]
Writes OUTPREFIX.a{l}.p{1,2}.ms and prints sizes.  A core row 'nonzero constant = 0' proves infeasibility directly."""
import sys, os, json, time, argparse
from fractions import Fraction
HERE = os.environ.get('NAG_RERUN_DIR', 'nag-rerun')
sys.path.insert(0, os.path.join(HERE, 'theirs'))
import numpy as np
from repair import build_from_recipe
from nagexport import sym_coefs
PR = [None]  # EXACT MODE (v-ex2-modp): one pass over Q

ap = argparse.ArgumentParser()
ap.add_argument('recipe'); ap.add_argument('out')
ap.add_argument('--keep', default=''); ap.add_argument('--mode', default='u')
ap.add_argument('--pins', default=''); ap.add_argument('--extra1', default=''); ap.add_argument('--extra2', default='')
ap.add_argument('--arrays', default='0'); ap.add_argument('--cap', type=int, default=30000)
ap.add_argument('--linsub', type=int, default=1)
ap.add_argument('--nonzero', default='')
ap.add_argument('--subs', default='')
ap.add_argument('--subsfile', default='')  # exact linsub: 'origidx=value,...'
args = ap.parse_args()

rec = json.load(open(args.recipe))
nsumm = [len(s) for (_, s, _) in rec['moves']]
orig_index = []                      # per move: list of original alpha indices
base = 0
for n in nsumm:
    orig_index.append(list(range(base, base + n))); base += n
NORIG = base
keep = {t: list(range(n)) for t, n in enumerate(nsumm)}
if args.keep:
    for part in args.keep.split(';'):
        t, ks = part.split(':'); keep[int(t)] = [int(k) for k in ks.split(',')]
sub_moves = []
sub2orig = []                        # sub-skeleton alpha index -> original alpha index
for t, (kind, summ, Dm) in enumerate(rec['moves']):
    sub_moves.append((kind, [summ[k] for k in keep[t]], Dm))
    sub2orig += [orig_index[t][k] for k in keep[t]]
dropped = sorted(set(range(NORIG)) - set(sub2orig))
pins = {}
for kv in [x for x in args.pins.split(',') if x]:
    k, v = kv.split('='); pins[int(k)] = Fraction(v)
tomat = lambda e: (np.array(e[0]) + 1j * np.array(e[1])) / 256.0
recipe = [(kind, [(src, tomat(b)) for (src, b) in summ], None if Dm is None else tomat(Dm)) for (kind, summ, Dm) in sub_moves]
t0 = time.time()
sk = build_from_recipe(rec['m'], rec['W'], rec['H'], recipe, 1)
img_sub = [len(mv['images']) for mv in sk.moves]
img_orig = [len(x) for x in rec['theirs']['images']] if 'theirs' in rec else None
cost = sum(1 + i for i in img_sub)
print(f"{rec['name']} keep={args.keep or 'all'} mode={args.mode} pins={args.pins or '-'} | gens {len(sk.gens)} images {img_sub}"
      f" (original {img_orig}) cost {cost} | dropped alphas {dropped}", flush=True)
if img_orig is not None and img_sub != img_orig:
    print("WARNING: image count differs from the original skeleton (generator indices of later moves would shift)", flush=True)
if args.mode == 'u':
    sk.fixed = []
first_of_move = set()
if args.mode == 'n':
    first_of_move = {sub2orig[a] for a in sk.fixed}
ok = sk.setup_final(args.cap)
print('closure', [(nr, nc) for (nr, nc, E) in sk.final] if ok else f'EXCEEDED cap {args.cap}', round(time.time() - t0, 1), 's', flush=True)
if not ok:
    sys.exit(3)
pc = sym_coefs(sk)                     # polynomials in SUB alpha indices (fixed ones already = 1 in mode n)

def read_extra(path, p):
    """alpha-only elements of an msolve GB file, as dicts {monomial(orig alpha idx...): coef mod p}"""
    if not path or not os.path.exists(path):
        return []
    txt = open(path).read()
    if 'TIMEOUT' in txt:
        return []
    body = txt[txt.index('['):].strip().rstrip(':').strip()[1:-1]
    out = []
    for poly in body.split(','):
        poly = poly.strip()
        if not poly or 'b' in poly:
            continue
        P = {}
        for term in poly.split('+'):
            fs = term.split('*'); c = int(fs[0]) % p; mono = []
            for f in fs[1:]:
                v, e = f.split('^') if '^' in f else (f, '1')
                mono += [int(v[1:]) - 1] * int(e)
            P[tuple(sorted(mono))] = (P.get(tuple(sorted(mono)), 0) + c) % p
        out.append(P)
    return out

def subst(P, p):
    """substitute dropped (=0), pinned and chart-fixed (=1) ORIGINAL alphas in a polynomial over orig indices"""
    out = {}
    for mono, c in P.items():
        c2 = c; m2 = []
        for a in mono:
            if a in dropped:
                c2 = 0; break
            if a in pins:
                v = pins[a]; c2 = c2 * v.numerator % p * pow(v.denominator, p - 2, p) % p
            elif a in first_of_move:
                pass
            else:
                m2.append(a)
        if c2:
            k = tuple(sorted(m2)); out[k] = (out.get(k, 0) + c2) % p
    return {k: v for k, v in out.items() if v}

def expo(k, names):
    """monomial as msolve factors WITH EXPONENTS: msolve 0.6.5 parses 'x*x' as x (found 14 Oct, nag-rerun), so a
    repeated variable must be written x^2"""
    out = []
    for v in sorted(set(k)):
        e = k.count(v); out.append(names[v] + (f"^{e}" if e > 1 else ""))
    return out
arrays = [int(x) for x in args.arrays.split(',')]
NZ = {int(x) for x in args.nonzero.split(',') if x}   # ORIGINAL alpha indices declared nonzero (t_k a_k = 1 added)
for jp, p in enumerate(PR):
    extras = []
    extras = [P for P in extras if P]
    if any(set(P) == {()} for P in extras):
        print(f'p{jp+1}: an extra (relaxation) equation became a nonzero constant in this chart -> chart infeasible')
    # --- linear substitution (sound: the linear alpha-only relaxation elements hold at every unit's alpha):
    # solve the degree-1 extras mod p for some alphas as affine forms in the others and substitute everywhere
    def pmulp(P, Q):
        out = {}
        for a, x in P.items():
            for b2, y in Q.items():
                k = tuple(sorted(a + b2)); out[k] = (out.get(k, 0) + x * y)
        return {k: v for k, v in out.items() if v}
    subs = {}
    def apply_subs(P):
        out = {}
        for mono, c in P.items():
            term = {(): c}
            for a in mono:
                term = pmulp(term, subs[a]) if a in subs else pmulp(term, {(a,): 1})
            for k, v in term.items():
                out[k] = (out.get(k, 0) + v)
        return {k: v for k, v in out.items() if v}
    for kv in [x for x in args.subs.split(',') if x]:
        kk, vv = kv.split('='); subs[int(kk)] = {(): Fraction(vv)}
    if args.subsfile:
        for kk, e in json.load(open(args.subsfile)).items():
            subs[int(kk)] = {tuple(int(u) for u in m.split(',') if u != ''): Fraction(c) for m, c in e.items()}
    if False:
        for P in [P for P in extras if all(len(k) <= 1 for k in P)]:
            P = apply_subs(P)
            vs = sorted({k[0] for k in P if k})
            if not vs:
                continue
            v = vs[-1]; inv = pow(P[(v,)], p - 2, p)
            expr = {k: (-x * inv) % p for k, x in P.items() if k != (v,)}
            expr = {k: x for k, x in expr.items() if x}
            for w in list(subs):
                subs[w] = apply_subs(subs[w]) if any(v in k for k in subs[w]) else subs[w]
            subs[v] = expr
            for w in list(subs):
                if w != v and any(v in k for k in subs[w]):
                    subs[w] = apply_subs(subs[w])
        extras = [apply_subs(P) for P in extras]
        extras = [P for P in extras if P]
        print(f'p{jp+1}: linear substitutions from the relaxation:', {f'a{v+1}': {('*'.join(f'a{u+1}' for u in k) or '1'): x for k, x in e.items()} for v, e in subs.items()}, flush=True)
        if any(set(P) == {()} for P in extras):
            print(f'p{jp+1}: an extra equation became a nonzero constant -> chart infeasible')
    def conv(Pc):
        out = {}
        for mono, co in Pc.items():
            assert abs(co.imag) < 1e-9 and abs(co.real - round(co.real)) < 1e-9
            c = Fraction(int(round(co.real))); m2 = []
            for a in mono:
                o = sub2orig[a]
                if o in pins:
                    raise SystemExit('pins not supported')
                else:
                    m2.append(o)
            if c:
                k = tuple(sorted(m2)); out[k] = (out.get(k, 0) + c)
        out = {k: v for k, v in out.items() if v}
        return apply_subs(out) if subs else out
    convc = {}
    for l in arrays:
        if l >= len(sk.final):
            continue
        nr, nc, E = sk.final[l]
        rows = [dict() for _ in range(nr)]
        for (r, c, i, j) in E:
            key = (i, j)
            if key not in convc:
                convc[key] = conv(pc[i][j])
            P = convc[key]
            R = rows[r]; Q = R.get(c, {})
            Q = dict(Q)
            for k, v in P.items():
                w = (Q.get(k, 0) + v)
                if w: Q[k] = w
                else: Q.pop(k, None)
            if Q: R[c] = Q
            else: R.pop(c, None)
        rows[0]['const'] = {(): Fraction(-1)}
        # constant-pivot elimination mod p (Markowitz-ish pivot choice)
        active = {i: R for i, R in enumerate(rows) if R}
        colrows = {}
        for i, R in active.items():
            for c in R:
                if c != 'const':
                    colrows.setdefault(c, set()).add(i)
        steps = 0
        def elig(R):
            return [c for c, P in R.items() if c != 'const' and len(P) == 1 and all(v in NZ for v in next(iter(P)))]
        fast = len(active) > 1500           # large systems: greedy pivot queue instead of a global Markowitz search
        queue = [i for i, R in active.items() if elig(R)] if fast else None
        while True:
            best = None
            if fast:
                while queue and best is None:
                    i = queue.pop()
                    if i not in active: continue
                    cs = elig(active[i])
                    if cs:
                        c = min(cs, key=lambda cc: len(colrows[cc])); best = (0, i, c)
                if best is None:
                    queue = [i for i, R in active.items() if elig(R)]     # rows that gained a pivot since
                    if not queue: break
                    continue
            else:
                for i, R in active.items():
                    for c, P in R.items():
                        if c != 'const' and len(P) == 1 and all(v in NZ for v in next(iter(P))):
                            score = (len(R) - 1) * (len(colrows[c]) - 1)
                            if best is None or score < best[0]:
                                best = (score, i, c)
                                if score == 0: break
                    if best is not None and best[0] == 0: break
            if best is None:
                break
            _, i, c = best
            Ri = active.pop(i)
            for cc in Ri:
                if cc != 'const': colrows[cc].discard(i)
            (k0, x0), = Ri[c].items()
            inv = Fraction(1) / x0
            for r2 in list(colrows[c]):
                R = active[r2]
                f = {k: (-v * inv) for k, v in R[c].items()}
                if k0:   # monomial pivot (nonzero alphas only): R := k0 * R - (R[c]/x0) * Ri, an equivalence when k0 != 0
                    R = {cc: {tuple(sorted(k + k0)): v for k, v in P.items()} for cc, P in R.items()}
                    active[r2] = R
                for cc, P in Ri.items():
                    Q = R.get(cc, {})
                    Q = dict(Q)
                    for k1, v1 in f.items():
                        for k2, v2 in P.items():
                            k = tuple(sorted(k1 + k2)); w = (Q.get(k, 0) + v1 * v2)
                            if w: Q[k] = w
                            else: Q.pop(k, None)
                    if Q:
                        if cc not in R and cc != 'const': colrows.setdefault(cc, set()).add(r2)
                        R[cc] = Q
                    else:
                        if cc in R and cc != 'const': colrows[cc].discard(r2)
                        R.pop(cc, None)
                assert c not in R
                if not R:
                    del active[r2]
            del colrows[c]
            steps += 1
        core = list(active.values())
        contra = [R for R in core if set(R) == {'const'} and set(R['const']) == {()}]
        cols = sorted({c for R in core for c in R if c != 'const'})
        avars = sorted({v for R in core for P in R.values() for k in P for v in k} | {v for P in extras for k in P for v in k})
        print(f"p{jp+1} array {l}: unit pivots {steps} | core rows {len(core)} beta cols {len(cols)} alphas {len(avars)}"
              + (" | CONTRADICTION: nonzero constant = 0 -> infeasible for EVERY alpha in this chart" if contra else "")
              + f" | extras {len(extras)}", flush=True)
        names = {v: f"a{v+1}" for v in avars}
        bn = {c: f"c{k+1}" for k, c in enumerate(cols)}
        lines = []
        for R in core:
            terms = []
            for c, P in R.items():
                for k, x in sorted(P.items()):
                    mon = expo(k, names) + ([] if c == 'const' else [bn[c]])
                    terms.append(f"{x}*{'*'.join(mon)}" if mon else f"{x}")
            if terms: lines.append("+".join(terms))
        for P in extras:
            terms = [f"{x}*{'*'.join(expo(k, names))}" if k else f"{x}" for k, x in sorted(P.items())]
            lines.append("+".join(terms))
        nzv = sorted(v for v in NZ if v in names)
        for v in nzv:
            lines.append(f"1*t{v+1}*{names[v]}+{p-1}")
        allv = [names[v] for v in avars] + [bn[c] for c in cols] + [f"t{v+1}" for v in nzv]
        if not lines:
            # every row eliminated by constant pivots: the system is solvable for EVERY alpha in this chart
            open(f"{args.out}.a{l}.Q.txt", 'w').write("# EMPTY CORE\n")

            print(f"p{jp+1} array {l}: EMPTY CORE -> solvable for every alpha in this chart", flush=True)
            continue
        open(f"{args.out}.a{l}.Q.txt", 'w').write(",".join(allv) + "\nQ\n" + ",\n".join(lines) + "\n")
print('done', round(time.time() - t0, 1), 's', flush=True)
