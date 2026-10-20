"""lcensus.py (w7-b5): exact F_2 census of the second-order fate of every producer of the class column's LEFTOVER
L_2 (array 2, frame A1, layer lam1 = 3, predecessor N_a) in the split-trunk 2 in 7 (and the 2-in-8: DEC8), pure Python.
Producers of L_2 (producer rule; t1/t2-injections, no side terms):
  Z1   z-column through a t1-stream term AT A1 (tau(M0) = Q, tau(N_c) = A1; layer 3): W3 Lemma L3's route.  The
       z-column b with b g = L_2 has main rows b h_c' pi_l at N_c' (pred Hbar_c'(W), W = Hbar_c^-1 tau^-1 (N_a)).
  Y2   y-column through a t2-stream companion at N_a (tau'(L0) = Q, tau'(N_c) = N_a; layer d(L0, N_c) = 2 =
       d(Q, N_a): needs profile (4,2)): t3 y-column of direction Hbar_c^-1 tau'^-1 (A1).
  Y1   y-column through a t1-stream companion at N_a (tau(N_c) = N_a; layer 3 = d(Q, N_a): profile (4,3)).
  Y12  y-column through a t1-stream's t2 term at N_a (tau(N1) = N_a; layer 2): t2 y-column of direction
       Hbar_2^-1 tau^-1 (A1), main rows at tau^-1(A1) for both arrays (pred N1).
Statuses as in e1census.py.  usage: lcensus.py seed nconf nsamp profA profB [8]"""
import sys, random, time
from collections import Counter
import e1census as C
from e1census import config, status, verdict, M0
from sp2 import L0, D as DS, mul, inv, image, dist
from fastP import fast_find_P as fp


def tau_t1(rng, cf, X, Y):
    """x2 t1-transport: tau(M0) = Q, tau(X) = Y"""
    return C.tau_with2(rng, cf, 1, X, Y)


def tau_t2(rng, cf, X, Y):
    """x2 t2-transport: tau'(L0) = Q, tau'(X) = Y.  tau' = qbar delta^-1 c2bar; with c2bar in P, write
    tau' = g h with g L0 = Q (g = qbar) and h in P: h(X) = g^-1 Y."""
    Q = cf['Q']
    if dist(L0, X) != dist(Q, Y):
        return None
    a2 = fp(rng, M0, image(DS, Q))
    g = mul(DS, mul(a2, DS))                      # g L0 = Q
    h = fp(rng, X, image(inv(g), Y))
    if h is None:
        return None
    t = mul(g, h)
    assert image(t, L0) == Q and image(t, X) == Y
    return t


def rows3(cf, X, only=None):
    out = []
    for c in ('a', 'b'):
        if C.DEC8 and only is not None and c != only:
            continue
        F = image(cf['H'][c], X); Nc = cf['N' + c]
        for l in (0, 1):
            out.append(status(cf, l, F, Nc))
    return out


def main():
    a = sys.argv[1:]
    seed, nconf, nsamp = int(a[0]), int(a[1]), int(a[2])
    pp = lambda s: None if s == '0' else (int(s[0]), int(s[1]))
    profA, profB = pp(a[3]), pp(a[4])
    C.DEC8 = len(a) > 5 and a[5] == '8'
    rng = random.Random(seed); t0 = time.time()
    out = Counter()
    for ci in range(nconf):
        cf = config(rng, profA, profB, 3, None)
        A1, Na = cf['A1'], cf['Na']
        for c in ('a', 'b'):
            Nc = cf['N' + c]
            # Z1: tau(N_c) = A1
            seen = set()
            for _ in range(nsamp):
                t = tau_t1(rng, cf, Nc, A1)
                if t is None:
                    break
                W = image(inv(cf['H'][c]), image(inv(t), Na))
                if W in seen:
                    continue
                seen.add(W)
                sv = []
                for cp in ('a', 'b'):
                    if C.DEC8 and cp != c:
                        continue
                    Ncp = cf['N' + cp]; Pr = image(cf['H'][cp], W)
                    for l in (0, 1):
                        anc = cf['P1'] if l == 0 else cf['Q']
                        lam = dist(anc, Ncp)
                        if lam <= 2:
                            sv.append('LOW')
                        elif lam == 4:
                            if Pr == cf['N1'] or (Pr in (cf['Na'], cf['Nb']) and dist(anc, Pr) >= 3):
                                sv.append('AUTO')
                            elif dist(anc, Pr) == 3:
                                sv.append('NEED')
                            else:
                                sv.append('DEAD')
                        else:
                            sv.append('NEED')
                out[('Z1 via ' + c, verdict(tuple(sv)))] += 1
            # Y2: t2-companion at N_a through branch c
            seen = set()
            for _ in range(nsamp):
                t = tau_t2(rng, cf, Nc, Na)
                if t is None:
                    break
                X = image(inv(cf['H'][c]), image(inv(t), A1))
                if X in seen:
                    continue
                seen.add(X)
                out[('Y2 via ' + c, verdict(tuple(rows3(cf, X, c))))] += 1
            # Y1: t1-companion at N_a through branch c
            seen = set()
            for _ in range(nsamp):
                t = tau_t1(rng, cf, Nc, Na)
                if t is None:
                    break
                X = image(inv(cf['H'][c]), image(inv(t), A1))
                if X in seen:
                    continue
                seen.add(X)
                out[('Y1 via ' + c, verdict(tuple(rows3(cf, X, c))))] += 1
        # Y12: t1-stream's t2 term at N_a
        seen = set()
        for _ in range(nsamp):
            t = tau_t1(rng, cf, cf['N1'], Na)
            if t is None:
                break
            R = image(inv(t), A1)
            if R in seen:
                continue
            seen.add(R)
            out[('Y12', verdict((status(cf, 0, R, cf['N1']), status(cf, 1, R, cf['N1']))))] += 1
        print("conf %d prof %s lam %s (%.0fs)" % (ci, cf['prof'], cf['lam'], time.time() - t0), flush=True)
    for k, v in sorted(out.items()):
        print("  ", k, v)


if __name__ == "__main__":
    main()
