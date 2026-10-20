"""kmax.py -- e3-bpexact: with kappa FREE (any combination of N letters of B'_m), how large can the top-alone kernel block
be in character a and irrep sigma?  For a generic row block w in Wsp (structured as in d4gen), w' = w R_sigma(e_a):
    Ssp = span over letters k of  w' R_sigma(k) |_(constrained slots + d0)
    Kach = { Z : [0 | Z] in Ssp }   (kernel blocks reachable with no junk on the constrained slots)
r(w) = max rank in Kach (generic element).  r(w) is a LOWER bound on r_a^max (a special w could do better).
The F2 criterion (NOTE s.4): F2 twisted by pi serves character a iff r_a >= d_sigma - m_a(pi); summed over a this needs
sum_a (d_sigma - r_a^max) <= d_sigma.  argv: seed psi N sigma(V|gl14|gl15|sq16|aff16) mode(all30|d4) [chars, default 1,0; r_1 = .. = r_6 by Galois]"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4gen import *
t0 = time.time()
seed, psi, N, sname, mode = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
L = Lift(seed=seed, psi=psi)
rng = random.Random(seed + 7); nrng = np.random.default_rng(seed + 9)
if sname == 'V': sig, ds = sig_nat, 16
else: sig, ds = make_sig(sname)
base = [b for (b, s) in L.Bm]
Bm = list(base); seen = set(key(b) for b in Bm)
# for representations trivial on the Paulis only the image modulo Paulis matters: dedupe on sigma(k) for those
seen_s = set()
mi = minv(L.m)
D = D4(L); P = L.P
ks = []
while len(ks) < N:
    w = mm(rng.choice(base), rng.choice(Bm)) if len(Bm) > len(base) or True else None
    if key(w) in seen: continue
    seen.add(key(w)); Bm.append(w)
    k = mm(L.m, w, mi)
    if sname in ('gl14', 'gl15'):          # these factor through the symplectic image
        sk = sigma(k)
        if sk in seen_s: continue
        seen_s.add(sk)
    ks.append(k)
acts = [D.act(k) for k in ks]
print(f'[kmax] sigma {sname} (dim {ds}), {len(ks)} letters of B\'_m, act computed {time.time() - t0:.0f}s', flush=True)
letters = [b for (b, s) in L.letters]
Bd1 = [b for b in letters if D.act(b)[D.d1][0] == D.d1]
rr = random.Random(5)
while len(Bd1) < 300: Bd1.append(mm(rr.choice(Bd1), rr.choice(Bd1)))
if sname == 'V':
    Alg = None; t = None
else:
    Alg = alg_basis([sig(D.act(b)[D.d1][1]) for b in Bd1], rr)
    t = [None] * 30
    for b in letters + Bd1:
        dd, beta = D.act(b)[D.d1]
        if t[dd] is None: t[dd] = sig(beta)
Fr_fixed = [d for d in range(30) if all(act[d][0] == d for act in acts[:300])]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(Fr_fixed))
cons = [d for d in range(30) if d != D.d0 and d not in Fr] if mode == 'all30' else D.Sx
blocks = cons + [D.d0]; pos = {d: i for i, d in enumerate(blocks)}
Rp = [D.R(p, sig, ds) for p in L.p]
out = {}
achars = [int(c) for c in (sys.argv[6] if len(sys.argv) > 6 else '1,0').split(',')]   # r_1 = ... = r_6 by Galois symmetry over Q(i)
for a in achars:
    if Alg is None:
        w = nrng.standard_normal((ds, 30 * ds)) + 1j * nrng.standard_normal((ds, 30 * ds))
    else:
        w = np.zeros((ds, 30 * ds), dtype=complex)
        for dd in range(30):
            Z = sum((nrng.standard_normal() + 1j * nrng.standard_normal()) * A for A in Alg)
            w[:, dd * ds:(dd + 1) * ds] = Z @ t[dd]
    Ea = sum(L.chars[a][j] * Rp[j] for j in range(7)) / 7
    wp = (w @ Ea).reshape(ds, 30, ds)
    S = np.zeros((len(ks), ds * len(blocks) * ds), dtype=complex)
    for j, act in enumerate(acts):
        o = np.zeros((ds, len(blocks), ds), dtype=complex)
        for d, (dd, beta) in enumerate(act):
            if dd in pos: o[:, pos[dd], :] += wp[:, d, :] @ sig(beta)
        S[j] = o.reshape(-1)
    nj = ds * len(cons) * ds
    S3 = S.reshape(len(ks), ds, len(blocks), ds)
    Sj = S3[:, :, :-1, :].reshape(len(ks), -1); Sk = S3[:, :, -1, :].reshape(len(ks), -1)
    # combinations c (rows) with c Sj = 0:
    U, s, Vh = np.linalg.svd(Sj, full_matrices=True)
    r = int((s > 1e-9 * s[0]).sum())
    Nc = U[:, r:].conj().T                      # c with c Sj = 0
    Kb = Nc @ Sk                                # reachable kernel blocks (as vectors)
    if Kb.shape[0] == 0:
        rk = 0
    else:
        z = (nrng.standard_normal(Kb.shape[0]) + 1j * nrng.standard_normal(Kb.shape[0])) @ Kb
        sv = np.linalg.svd(z.reshape(ds, ds), compute_uv=False)
        rk = int((sv > 1e-8 * np.abs(S).max()).sum())
    Ep = sum(L.chars[a][j] * sig(L.p[j]) for j in range(7)) / 7
    ma = int(round(np.trace(Ep).real))
    spanS = int((np.linalg.svd(S, compute_uv=False) > 1e-9 * s[0]).sum())
    out[a] = (rk, ma)
    print(f'[kmax] {sname} {mode} a={a}: span of w\'R(k) {spanS} (of {S.shape[1]}), junk-free combos {Nc.shape[0]}, '
          f'reachable kernel-block rank r = {rk} (d_sigma {ds}); m_a(p) = {ma}; F2(p) criterion r >= {ds - ma}: '
          f'{"ok" if rk >= ds - ma else "FAIL"}  {time.time() - t0:.0f}s', flush=True)
deficit = sum(ds - out[a][0] for a in achars) + (5 * (ds - out[1][0]) if 1 in out and len(achars) == 2 else 0)
print(f'[kmax] {sname} {mode}: total deficit sum_a (d - r_a) = {deficit}  vs  d_sigma = {ds}  (any F2 twist needs deficit <= d)')
print('done %.0fs' % (time.time() - t0))
