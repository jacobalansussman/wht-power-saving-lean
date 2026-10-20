#!/usr/bin/env python3
"""an2-smallnet / designs.py (own code, stdlib + fractions; python3 -I).  For the port designs a small-block crossing
network could use (root systems A_n, D_n, E6, E7, E8, the hub's P112, written from scratch), compute exactly:
  n ports (lines), t = rank over Q of the line Gram A (the centres an index code needs at least), the binary label
  Gram B = J + G mod 2 (junk only on label-orthogonal pairs), the window k = rank_F2 B, omega = the largest set of
  mutually label-non-orthogonal ports (= mutually orthogonal lines; the floor t >= omega), junk per port s, and
  sigma = min over ports of rank_F2 of the labels of its junk partners: the frame levels every target must take up
  in every call (it meets each junk partner's source at that source's frame, which holds the partner's label);
  by time-reversal duality every source takes up as many.  Block floor of a call layout in which some array takes
  part in 3 calls (shears.py: unavoidable): m >= 1 + 3 sigma.  Hub accounting copy ratio t (k - 1) / n."""
import itertools, json, os, sys
from fractions import Fraction

def roots_E8():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * 8; v[i] = 2 * si; v[j] = 2 * sj; R.append(tuple(v))   # doubled coordinates
    for signs in itertools.product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            R.append(tuple(signs))
    return R   # norm 8 in doubled coordinates

def ip(a, b): return sum(x * y for x, y in zip(a, b))

def lines(R):
    seen, L = set(), []
    for r in R:
        key = max(r, tuple(-x for x in r))
        if key not in seen:
            seen.add(key); L.append(key)
    return L

def roots_D(n):
    R = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [0] * n; v[i] = si; v[j] = sj; R.append(tuple(v))
    return R

def roots_A(n):
    R = []
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                v = [0] * (n + 1); v[i] = 1; v[j] = -1; R.append(tuple(v))
    return R

def rank_Q(rows):
    M = [[Fraction(x) for x in r] for r in rows]
    rk, col, ncol = 0, 0, len(M[0]) if M else 0
    for col in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][col] != 0:
                f = M[i][col] / M[rk][col]
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk

def rank_F2(rows):
    basis = []
    for r in rows:
        v = int(''.join('1' if x % 2 else '0' for x in r), 2) if isinstance(r, (list, tuple)) else r
        for b in basis:
            v = min(v, v ^ b)
        if v: basis.append(v)
    return len(basis)

def max_clique(adj, n):
    best = []
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    def expand(R, P):
        nonlocal best
        if len(R) + len(P) <= len(best): return
        if not P:
            if len(R) > len(best): best = R[:]
            return
        for v in list(P):
            if len(R) + len(P) <= len(best): return
            expand(R + [v], P & adj[v])
            P = P - {v}
    expand([], set(order))
    return best

def analyse(name, L):
    n = len(L)
    norm = ip(L[0], L[0])
    G = [[ip(a, b) for b in L] for a in L]            # integer Gram
    t = rank_Q(G)
    # label Gram over F_2: B = 1 + <r,s>_norm2  (orthogonal or equal -> 1, junk -> 0); <r,s> in units of norm/2
    unit = norm // 2
    B = [[(1 + (G[i][j] // unit)) % 2 if G[i][j] % unit == 0 else None for j in range(n)] for i in range(n)]
    assert all(B[i][j] is not None for i in range(n) for j in range(n)), 'not a two-angle design'
    k = rank_F2(B)
    junk = [[j for j in range(n) if j != i and G[i][j] != 0] for i in range(n)]
    s = len(junk[0])
    sig = min(rank_F2([B[j] for j in junk[i]]) for i in range(n))
    sig_max = max(rank_F2([B[j] for j in junk[i]]) for i in range(n))
    adj = [set(j for j in range(n) if j != i and B[i][j] == 1) for i in range(n)]
    om = len(max_clique(adj, n))
    tight = all(sum(Fraction(G[i][l] * G[l][j], norm * norm) for l in range(n)) == Fraction(n, t) * Fraction(G[i][j], norm)
                for i in range(0, n, max(1, n // 7)) for j in range(n))
    cr = Fraction(t * (k - 1), n)
    rec = dict(design=name, n=n, t=t, omega=omega_ok(om, t), k=k, junk_per_port=s, sigma=sig, sigma_max=sig_max,
               tight=tight, n_over_t=float(Fraction(n, t)), copy_ratio_hub=float(cr),
               block_floor_3calls=1 + 3 * sig, block_floor_2calls=1 + 2 * sig)
    print(f"{name:6s} n={n:4d} t={t:2d} omega={om:2d} window k={k:2d} junk/port={s:3d} sigma={sig}..{sig_max} tight={tight}"
          f" n/t={float(Fraction(n,t)):6.3f} copy ratio t(k-1)/n={float(cr):.3f}  block floor: 1+3sigma={1+3*sig}, 1+2sigma={1+2*sig}",
          flush=True)
    return rec

def omega_ok(om, t):
    return om

def main():
    out = []
    E8 = lines(roots_E8())
    frame = set()
    for i in range(0, 8, 2):
        for sj in (1, -1):
            v = [0] * 8; v[i] = 2; v[i + 1] = 2 * sj; frame.add(tuple(v))
    P112 = [l for l in E8 if l not in frame]
    r0 = (2, 2, 0, 0, 0, 0, 0, 0)
    E7 = [l for l in E8 if ip(l, r0) == 0]
    r1 = (2, -2, 0, 0, 0, 0, 0, 0); r2 = (0, 2, -2, 0, 0, 0, 0, 0)          # an A2 inside E8 (two roots at 120 deg)
    r2 = (0, 2, 2, 0, 0, 0, 0, 0)
    E6 = [l for l in E8 if ip(l, (2, -2, 0, 0, 0, 0, 0, 0)) == 0 and ip(l, (0, 2, -2, 0, 0, 0, 0, 0)) == 0]
    designs = [('E8', E8), ('P112', P112), ('E7', E7), ('E6', E6)]
    for n_ in (4, 5, 6, 7, 8):
        designs.append((f'D{n_}', lines(roots_D(n_))))
    for n_ in (3, 4, 5, 6, 7):
        designs.append((f'A{n_}', lines(roots_A(n_))))
    for name, L in designs:
        out.append(analyse(name, L))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'designs.json'), 'w'), indent=1)

if __name__ == '__main__':
    main()
