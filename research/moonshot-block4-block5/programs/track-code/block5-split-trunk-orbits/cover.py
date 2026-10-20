"""cover.py (w7-b5): two exact F_2 checks for NOTE.md sections 2-3, pure Python.

(A) Lemma E1-sigma, exhaustive form: for k' in {a, b}, the set {sigma(A_o) : sigma in Stab(M0), sigma(N_own) = N_k'}
    equals the set of ALL neighbours Y of N_k' with d(M0, Y) = 4 (transitivity check by sampling until the sampled
    image set stops growing and comparing with the full set), and the second-order status of the E1 producer for EVERY
    such Y (2-in-7 split trunk: rows on both branches; 2-in-8: rows on branch k' only).
(B) Lemma T (trunk cover): for a t1-injection of x_i (transport tau, tau(M0) = anchor_i), its two own-branch Lemma-Q rows
    R1 (output 1's class column: frame tau(A1), pred tau(N_a)) and R2 (output 2's: tau(A2), pred tau(N_b)) are both at
    layer 4.  Force each single trunk escape exactly (E3: tau(N) = N1; E2: tau(N) = N_c with d(anchor_i, N_c) = 3;
    E4: tau(A) = N_c with d(anchor_i, N_c) = 4) for one row by sampling the coset {tau : tau(M0) = anchor, tau(X) = Y},
    and record how often the OTHER row also has a trunk escape.
usage: cover.py seed nconf nsamp [profA profB lam1 lam2]"""
import sys, random, time
from collections import Counter
import e1census as C
from e1census import config, status, verdict, tau_with2, M0
from sp2 import L0, D as DS, mul, inv, image, dist, neighbors
from fastP import fast_find_P as fp


def rows_status(cf, X, only):
    out = []
    for c in ('a', 'b'):
        if only is not None and c != only:
            continue
        F = image(cf['H'][c], X); Nc = cf['N' + c]
        for l in (0, 1):
            out.append(status(cf, l, F, Nc))
    return tuple(out)


def partA(rng, cf, nsamp, out):
    for o, bo in ((0, 'a'), (1, 'b')):
        Ao = cf['A1'] if o == 0 else cf['A2']
        Zo = image(inv(cf['H'][bo]), Ao)
        for kp in ('a', 'b'):
            Nk = cf['N' + kp]
            full = {Y for Y in neighbors(Nk) if dist(M0, Y) == 4}
            seen = set()
            Xa = image(DS, cf['N' + bo]); Xk = image(DS, Nk)
            for _ in range(nsamp):
                u = fp(rng, Xa, Xk)
                sig = mul(DS, mul(u, DS))
                seen.add(image(sig, Ao))
            out['transitive'][(len(seen), len(full), seen == full)] += 1
            for Y in full:
                X = image(inv(cf['H'][kp]), Y)
                tw = '/twin' if X == Zo else ''
                out['E1-7'][verdict(rows_status(cf, X, None)) + tw] += 1
                out['E1-8'][verdict(rows_status(cf, X, kp)) + tw] += 1


def escapes(cf, i, tau, which):
    """trunk escapes of row R1 (which=0) or R2 (which=1) of the x_i stream with transport tau"""
    anc = cf['P1'] if i == 0 else cf['Q']
    A = cf['A1'] if which == 0 else cf['A2']
    N = cf['Na'] if which == 0 else cf['Nb']
    F = image(tau, A); Pr = image(tau, N)
    assert dist(anc, F) == 4 and dist(anc, Pr) == 3
    out = []
    if Pr == cf['N1']:
        out.append('E3')
    for c in ('a', 'b'):
        Nc = cf['N' + c]
        if Pr == Nc and dist(anc, Nc) == 3:
            out.append('E2' + c)
        if F == Nc and dist(anc, Nc) == 4:
            out.append('E4' + c)
    return out


def partB(rng, cf, nsamp, out):
    for i in (0, 1):
        anc = cf['P1'] if i == 0 else cf['Q']
        for which in (0, 1):
            A = cf['A1'] if which == 0 else cf['A2']
            N = cf['Na'] if which == 0 else cf['Nb']
            forced = [('E3', N, cf['N1'])]
            for c in ('a', 'b'):
                Nc = cf['N' + c]
                if dist(anc, Nc) == 3:
                    forced.append(('E2' + c, N, Nc))
                if dist(anc, Nc) == 4:
                    forced.append(('E4' + c, A, Nc))
            for (name, X, Y) in forced:
                hit = 0; tot = 0
                for _ in range(nsamp):
                    t = tau_with2(rng, cf, i, X, Y)
                    if t is None:
                        break
                    tot += 1
                    if escapes(cf, i, t, 1 - which):
                        hit += 1
                out['T'][('x%d' % (i + 1), 'R%d by %s' % (which + 1, name), 'other row covered', hit, 'of', tot)] += 1


def main():
    a = sys.argv[1:]
    seed, nconf, nsamp = int(a[0]), int(a[1]), int(a[2])
    pp = lambda s: None if s == '0' else (int(s[0]), int(s[1]))
    profA = pp(a[3]) if len(a) > 3 else None
    profB = pp(a[4]) if len(a) > 4 else None
    lam1 = (None if a[5] == '0' else int(a[5])) if len(a) > 5 else None
    lam2 = (None if a[6] == '0' else int(a[6])) if len(a) > 6 else None
    rng = random.Random(seed); t0 = time.time()
    out = {'transitive': Counter(), 'E1-7': Counter(), 'E1-8': Counter(), 'T': Counter()}
    for ci in range(nconf):
        cf = config(rng, profA, profB, lam1, lam2)
        partA(rng, cf, nsamp, out)
        partB(rng, cf, nsamp, out)
        print("conf %d prof %s lam %s done (%.0fs)" % (ci, cf['prof'], cf['lam'], time.time() - t0), flush=True)
    for k, Cn in out.items():
        print(k)
        for kk, v in sorted(Cn.items(), key=str):
            print("   ", kk, v)


if __name__ == "__main__":
    main()
