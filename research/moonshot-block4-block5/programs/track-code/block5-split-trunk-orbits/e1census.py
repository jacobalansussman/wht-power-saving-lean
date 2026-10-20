"""e1census.py (w7-b5): exact F_2 census of the SECOND-ORDER fate of every escape route of a Lemma-Q pollution row
(b5-twobranch PROOF sec. 2b) in the split-trunk 2 in 7 at block 5 (engine game), pure Python.

Global picture (w3-b5 PROOF sec. 1): trunk frames M0 (s), L0 (t1), N1 (t2), branch frames N_c = Hbar_c(L0) (t3), anchors
P1, Q; targets A1 ~ N_a, A2 ~ N_b.  A y-column on t3 with DIRECTION X (a neighbour of L0) sends a t3 term with
theta = Hbar_c (main stream, branch c) to the row frame Hbar_c(X), predecessor N_c; a stream with transport tau to
tau Hbar_c(X), predecessor tau(N_c).

Output 1's class column (direction Z = Hbar_a^-1 A1) applied to an x2 t1-stream j (transport tau_j, tau_j M0 = Q) has the
Lemma-Q row r at tau_j(A1) (layer 4), predecessor tau_j(N_a) (layer 3).  Its possible extra producers:
  E1  stream k via branch k' with tau_k(N_k') = tau_j(N_a).  Lemma E1-sigma: sigma = tau_k^-1 tau_j is an ARBITRARY
      element of Stab(M0) with sigma(N_a) = N_k'; the producer is the t3 y-column with direction Hbar_k'^-1 sigma(A1);
      its main rows sit at Hbar_c Hbar_k'^-1 sigma(A1) (pred N_c), for every branch c and both arrays.  Independent of tau_j.
  E2  a branch main term h_c' q at tau_j(N_a) (needs d(Q, N_c') = 3): t3 y-column, direction Hbar_c'^-1 tau_j(A1).
  E3  t2's main term at tau_j(N_a) = N1: t2 y-column, main rows at tau_j(A1) for both arrays (pred N1).
  E4  a branch main term at tau_j(A1) (z-column on t3): main rows b'' h_c pi_l at N_c, pred Hbar_c(W),
      W = Hbar_c'^-1 tau_j(N_a).
For each producer, every MAIN row (both arrays, every branch) gets a status:
  FREE (on the array's own K-frame A_l), LOW (layer <= 2: irrelevant for the d >= 3 tests), AUTO (an automatic second
  producer from the trunk: z via t2 main at N1, z via a branch main term at the row's frame with the right layer),
  NEED (a second producer needs a further injected stream or companion), DEAD (no possible second producer at all: the
  row forces the producer to 0).  Rules (proved in NOTE.md sec. 2): for array l, row frame F, predecessor frame Pr,
  lam = d(anchor_l, F), nu = d(anchor_l, Pr):  lam = 5: FREE if F = A_l else DEAD.  lam = 4: AUTO if F is a branch
  frame with d(anchor_l, F) = 4; else NEED if nu = 3 (a companion stream at Pr) else DEAD.  lam = 3: AUTO if F = N1 or a
  branch frame with d(anchor_l, F) = 3; else NEED.  lam <= 2: LOW.
usage: e1census.py seed nconf nsig [profA profB lam1 lam2]"""
import sys, random, time
from collections import Counter
from sp2 import M, L0, D as DS, mul, inv, image, dist, neighbors
from fastP import fast_find_P as fp

M0 = image(DS, L0)
NLET = 12
FILTER = True
DEC8 = False


def config(rng, profA, profB, lam1, lam2):
    N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
    for _ in range(400):
        N1 = rng.choice(N1s)
        brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
        anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
        for _ in range(40):
            P1, Q = rng.sample(anchors, 2)
            if dist(P1, Q) != 2:
                continue
            pr = {N: (dist(P1, N), dist(Q, N)) for N in brs}
            Nas = [N for N in brs if pr[N][0] == 4 and (profA is None or pr[N] == profA)]
            Nbs = [N for N in brs if pr[N][1] == 4 and (profB is None or pr[N] == profB)]
            rng.shuffle(Nas); rng.shuffle(Nbs)
            for Na in Nas[:3]:
                A1s = [A for A in sorted(neighbors(Na)) if dist(P1, A) == 5 and dist(Q, A) <= 4 and (lam1 is None or dist(Q, A) == lam1)]
                if not A1s:
                    continue
                for Nb in Nbs[:3]:
                    if Nb == Na:
                        continue
                    A2s = [A for A in sorted(neighbors(Nb)) if dist(Q, A) == 5 and dist(P1, A) <= 4 and (lam2 is None or dist(P1, A) == lam2)]
                    if not A2s:
                        continue
                    A1 = rng.choice(A1s); A2 = rng.choice(A2s)
                    c2 = fp(rng, N1, M0)
                    for _l in range(NLET):
                        c3a = fp(rng, image(DS, image(c2, Na)), M0)
                        c3b = fp(rng, image(DS, image(c2, Nb)), M0)
                        Ha = mul(inv(c2), mul(DS, mul(inv(c3a), DS)))
                        Hb = mul(inv(c2), mul(DS, mul(inv(c3b), DS)))
                        assert image(Ha, L0) == Na and image(Hb, L0) == Nb and image(Ha, M0) == N1 and image(Hb, M0) == N1
                        H2 = mul(inv(c2), DS)
                        assert image(H2, L0) == N1
                        cf = dict(N1=N1, P1=P1, Q=Q, Na=Na, Nb=Nb, A1=A1, A2=A2, H={'a': Ha, 'b': Hb}, H2=H2,
                                  prof=(pr[Na], pr[Nb]), lam=(dist(Q, A1), dist(P1, A2)), c2=c2, c3a=c3a, c3b=c3b)
                        if not FILTER or not any(verdict(sv) == 'DEAD' for sv in class_ok(cf)):
                            return cf
    return None


def status(cf, l, F, Pr):
    anc = cf['P1'] if l == 0 else cf['Q']
    A = cf['A1'] if l == 0 else cf['A2']
    lam = dist(anc, F)
    brf = {cf['Na'], cf['Nb']}
    if lam == 5:
        return 'FREE' if F == A else 'DEAD'
    if lam == 4:
        if F in brf:
            return 'AUTO'
        return 'NEED' if dist(anc, Pr) == 3 else 'DEAD'
    if lam == 3:
        if F == cf['N1'] or F in brf:
            return 'AUTO'
        return 'NEED'
    return 'LOW'


def t3_rows(cf, X, only=None):
    """status vector of a t3 y-column with direction X: (array, branch) -> status.  DEC8 (the decoupled 2-in-8 shape
    s1 s2 t1 t2 ua ub e1 e2): a y-column on u_c has main rows on branch c only (`only`); the other branch's entries are
    reported as 'LOW' (absent)."""
    out = []
    for c in ('a', 'b'):
        if DEC8 and only is not None and c != only:
            out += ['LOW', 'LOW']; continue
        F = image(cf['H'][c], X); Nc = cf['N' + c]
        for l in (0, 1):
            out.append(status(cf, l, F, Nc))
    return tuple(out)


def verdict(sv):
    if 'DEAD' in sv:
        return 'DEAD'
    n = sv.count('NEED')
    return 'OK0' if n == 0 else 'NEED%d' % n


def tau_with(rng, cf, X, Y):
    """random tau with tau(M0) = Q and tau(X) = Y (None if impossible)"""
    Q = cf['Q']
    if dist(M0, X) != dist(Q, Y):
        return None
    a2 = fp(rng, M0, image(DS, Q))
    qb = mul(DS, mul(a2, DS))
    b = fp(rng, image(inv(qb), Y), image(DS, X))
    if b is None:
        return None
    t = mul(qb, mul(inv(b), DS))
    assert image(t, M0) == Q and image(t, X) == Y
    return t



def class_ok(cf):
    """the class columns' own main rows (other-branch images and leftovers) must not be DEAD"""
    out = []
    for o, bo in ((0, 'a'), (1, 'b')):
        Ao = cf['A1'] if o == 0 else cf['A2']
        Z = image(inv(cf['H'][bo]), Ao)
        sv = t3_rows(cf, Z, only=bo)
        out.append(sv)
    return out


def producers(rng, cf, o, i, nsig):
    """census of the second-order status of every escape route of the pollution row of output o's class column on an
    input-i t1-stream.  Returns dict route -> Counter(verdict)."""
    bo = 'a' if o == 0 else 'b'
    Ao = cf['A1'] if o == 0 else cf['A2']
    Nbo = cf['N' + bo]
    anc_i = cf['P1'] if i == 0 else cf['Q']
    Zo = image(inv(cf['H'][bo]), Ao)
    res = {}
    # E1
    for kp in ('a', 'b'):
        C = Counter(); seen = set()
        Xa = image(DS, Nbo); Xk = image(DS, cf['N' + kp])
        for _ in range(nsig):
            u = fp(rng, Xa, Xk)
            if u is None:
                break
            sig = mul(DS, mul(u, DS))
            Y = image(sig, Ao)
            if Y in seen:
                continue
            seen.add(Y)
            X = image(inv(cf['H'][kp]), Y)
            sv = t3_rows(cf, X, only=kp)
            C[verdict(sv) + ('/twin' if X == Zo else '')] += 1
        res['E1' + kp] = C
    # E2
    C = Counter()
    for cp in ('a', 'b'):
        Ncp = cf['N' + cp]
        if dist(anc_i, Ncp) != 3:
            continue
        seen = set()
        for _ in range(nsig):
            t = tau_with2(rng, cf, i, Nbo, Ncp)
            if t is None:
                break
            R = image(t, Ao)
            if R in seen:
                continue
            seen.add(R)
            X = image(inv(cf['H'][cp]), R)
            sv = list(t3_rows(cf, X, only=cp))
            sv[('a', 'b').index(cp) * 2 + i] = 'SELF'
            C[cp + ':' + verdict(tuple(sv))] += 1
    res['E2'] = C
    # E3
    C = Counter(); seen = set()
    for _ in range(nsig):
        t = tau_with2(rng, cf, i, Nbo, cf['N1'])
        if t is None:
            break
        R = image(t, Ao)
        if R in seen:
            continue
        seen.add(R)
        C[verdict((status(cf, 1 - i, R, cf['N1']),))] += 1
    res['E3'] = C
    # E4
    C = Counter()
    for cp in ('a', 'b'):
        Ncp = cf['N' + cp]
        if dist(anc_i, Ncp) != 4:
            continue
        seen = set()
        for _ in range(nsig):
            t = tau_with2(rng, cf, i, Ao, Ncp)
            if t is None:
                break
            W = image(inv(cf['H'][cp]), image(t, Nbo))
            if W in seen:
                continue
            seen.add(W)
            sv = []
            for c in ('a', 'b'):
                Nc = cf['N' + c]; Pr = image(cf['H'][c], W)
                if DEC8 and c != cp:
                    continue
                for l in (0, 1):
                    if c == cp and l == i:
                        continue
                    anc = cf['P1'] if l == 0 else cf['Q']
                    lam = dist(anc, Nc)
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
            C[cp + ':' + verdict(tuple(sv))] += 1
    res['E4'] = C
    return res


def tau_with2(rng, cf, i, X, Y):
    anc = cf['P1'] if i == 0 else cf['Q']
    if dist(M0, X) != dist(anc, Y):
        return None
    a2 = fp(rng, M0, image(DS, anc))
    qb = mul(DS, mul(a2, DS))
    b = fp(rng, image(inv(qb), Y), image(DS, X))
    if b is None:
        return None
    t = mul(qb, mul(inv(b), DS))
    assert image(t, M0) == anc and image(t, X) == Y
    return t


def alive(C):
    return any(not k.split(':')[-1].startswith('DEAD') for k in C)


def main():
    a = sys.argv[1:]
    seed, nconf, nsig = int(a[0]), int(a[1]), int(a[2])
    def pp(s):
        return None if s == '0' else (int(s[0]), int(s[1]))
    profA = pp(a[3]) if len(a) > 3 else None
    profB = pp(a[4]) if len(a) > 4 else None
    lam1 = (None if a[5] == '0' else int(a[5])) if len(a) > 5 else None
    lam2 = (None if a[6] == '0' else int(a[6])) if len(a) > 6 else None
    global DEC8
    DEC8 = len(a) > 7 and a[7] == '8'
    print("mode:", "decoupled 2-in-8 (ua, ub separate moves)" if DEC8 else "split-trunk 2 in 7", flush=True)
    rng = random.Random(seed); t0 = time.time()
    tot = {}; sets = Counter(); done = 0; skipped = 0
    while done < nconf:
        cf = config(rng, profA, profB, lam1, lam2)
        if cf is None:
            print("no config"); return
        co = class_ok(cf)
        if any(verdict(sv) == 'DEAD' for sv in co):
            skipped += 1
            if skipped > 50 * nconf:
                print("too many dead class columns"); break
            continue
        done += 1
        line = []
        for (o, i) in ((0, 1), (1, 0), (0, 0), (1, 1)):
            res = producers(rng, cf, o, i, nsig)
            key = 'out%d/x%d' % (o + 1, i + 1)
            al = tuple(sorted(r for r, C in res.items() if alive(C)))
            sets[(key, al)] += 1
            for r, C in res.items():
                tot.setdefault((key, r), Counter()).update(C)
            line.append("%s %s" % (key, {r: dict(C) for r, C in res.items() if C}))
        print("conf %d prof %s lam %s class rows %s (%.0fs)\n   " % (done, cf['prof'], cf['lam'], co, time.time() - t0)
              + "\n   ".join(line), flush=True)
    print("skipped configs with a DEAD class column: %d" % skipped)
    print("TOTALS (distinct producer images, summed over configs):")
    for k, C in sorted(tot.items()):
        print("  %s %s: %s" % (k[0], k[1], dict(sorted(C.items()))))
    print("(pollution kind, routes with a non-DEAD option): number of configs")
    for k, v in sorted(sets.items()):
        print("  %s: %d" % (k, v))


if __name__ == "__main__":
    main()
