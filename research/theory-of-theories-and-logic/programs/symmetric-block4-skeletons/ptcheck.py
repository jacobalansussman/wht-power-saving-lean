"""ptcheck.py: for each hit with an exact rational relaxation point in ex2-modp/out/relax_hits.txt:
 (1) find the alpha-index convention under which MY relaxation is solvable at that alpha (exact over Q, rank test),
 (2) at that alpha, test MY full closure system of array 0 (summands with alpha = 0 dropped), exact over Q via two primes.
A point where the full system were solvable would be a unit candidate."""
import sys, os, re, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import vnr
REC = os.path.join(HERE, 'recipes')
PS = [1073741789, 1073741741]

def fmod(x, p):
    x = Fraction(x); return x.numerator % p * pow(x.denominator % p, p - 2, p) % p

def relax_feasible(sk, al, p):
    rel = sk.relaxation()
    ng = len(sk.gens)
    M = np.zeros((len(rel), ng + 1), dtype=np.int64)
    alp = {a: fmod(v, p) for a, v in al.items()}
    for i, (R, c) in enumerate(rel):
        for s, P in R.items():
            M[i, s] = vnr.peval(P, alp, p)
        M[i, ng] = (-c) % p
    return vnr.rank_mod(M[:, :ng], p), vnr.rank_mod(M, p)

def full_feasible(rec, al, p):
    zero = {a for a, v in al.items() if v == 0}
    sk = vnr.Skel(rec, zero=zero)
    cl = sk.closure(0, 3000)
    if cl is None: return None
    A, nc = cl
    alp = {a: fmod(v, p) for a, v in al.items()}
    return (len(A), nc) + vnr.full_screen(A, nc, alp, p)

lines = open(sys.argv[1]).read().splitlines()
for i, l in enumerate(lines):
    m = re.match(r'(\S+)\.u: ', l)
    if not m: continue
    name = m.group(1)
    if not os.path.exists(os.path.join(REC, name + '.json')): continue
    nxt = lines[i + 1]
    if not nxt.startswith('EXISTS'):
        print(name, 'no exact point in ex2-modp'); continue
    asg = dict(re.findall(r'\b(a\d+)=(-?\d+(?:/\d+)?)', nxt.split(':', 1)[1]))
    rec = json.load(open(os.path.join(REC, name + '.json')))
    sk = vnr.Skel(rec)
    for off in (0, 1):
        al = {a: Fraction(0) for a in range(sk.na)}
        ok = True
        for k, v in asg.items():
            j = int(k[1:]) - off
            if j < 0 or j >= sk.na: ok = False; break
            al[j] = Fraction(v)
        if not ok: continue
        rk = [relax_feasible(sk, al, p) for p in PS]
        if all(r[0] == r[1] for r in rk):
            ff = [full_feasible(rec, al, p) for p in PS]
            verdict = 'FULL SOLVABLE (unit candidate!)' if all(f and f[2] == f[3] for f in ff) else 'full infeasible'
            print(name, 'offset', off, 'relaxation solvable at the point (both primes)', 'full closure', ff[0][:2],
                  'ranks', [f[2:] for f in ff], verdict, flush=True)
            break
    else:
        print(name, 'point does NOT satisfy my relaxation under either index convention', flush=True)
