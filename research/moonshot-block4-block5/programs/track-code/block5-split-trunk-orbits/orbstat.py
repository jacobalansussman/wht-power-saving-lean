"""orbstat.py (w7-b5): exact F_2 statistic for frame-closed sigma-orbits (sigma in Stab(M0) mapping {N_a, N_b} to itself:
'fix' fixes both, 'swap' swaps them).  For each own-branch Lemma-Q row type (output o on its class branch c_o) and each
j = 1..ord-1, the E1 producer inside the orbit has direction Hbar_k'^-1 sigma^j (A_o), k' = the branch sigma^j(N_c_o)
(Lemma E1-sigma); classify it TWIN (class direction: Corollary O), DEAD, or ALIVE (non-class, no DEAD main row).
A frame-closed orbit can escape Corollary O only if every own-branch row type has an ALIVE j.
usage: orbstat.py seed nconf ntries [profA profB lam1 lam2 [8]]"""
import sys, random, time
from collections import Counter
import e1census as C
from e1census import config, verdict, M0
from sp2 import D as DS, mul, inv, image, order as sorder
from fastP import fast_find_P as fp


def main():
    a = sys.argv[1:]
    seed, nconf, ntries = int(a[0]), int(a[1]), int(a[2])
    pp = lambda s: None if s == '0' else (int(s[0]), int(s[1]))
    profA = pp(a[3]) if len(a) > 3 else None
    profB = pp(a[4]) if len(a) > 4 else None
    lam1 = (None if a[5] == '0' else int(a[5])) if len(a) > 5 else None
    lam2 = (None if a[6] == '0' else int(a[6])) if len(a) > 6 else None
    C.DEC8 = len(a) > 7 and a[7] == '8'
    rng = random.Random(seed); t0 = time.time()
    tot = Counter(); per_sigma = Counter()
    for ci in range(nconf):
        cf = config(rng, profA, profB, lam1, lam2)
        Na, Nb = cf['Na'], cf['Nb']
        Xa, Xb = image(DS, Na), image(DS, Nb)
        Z = {0: image(inv(cf['H']['a']), cf['A1']), 1: image(inv(cf['H']['b']), cf['A2'])}
        for kind in ('fix', 'swap'):
            n = 0
            for t in range(ntries):
                u = fp(rng, Xa, Xa if kind == 'fix' else Xb)
                if image(u, Xb) != (Xb if kind == 'fix' else Xa):
                    continue
                o = sorder(u)
                if o is None or o > 60:
                    continue
                n += 1
                sig = mul(DS, mul(u, DS))
                ok_all = True; cls = {}
                for o_, co in ((0, 'a'), (1, 'b')):
                    Ao = cf['A1'] if o_ == 0 else cf['A2']
                    sj = sig; has_alive = False
                    for j in range(1, o):
                        Nk = image(sj, cf['N' + co]); kp = 'a' if Nk == Na else 'b'
                        X = image(inv(cf['H'][kp]), image(sj, Ao))
                        if X == Z[o_] and (not C.DEC8 or kp == co):
                            s = 'TWIN'
                        else:
                            s = 'DEAD' if verdict(C.t3_rows(cf, X, only=kp)) == 'DEAD' else 'ALIVE'
                        tot[(kind, 'out%d' % (o_ + 1), s)] += 1
                        has_alive |= (s == 'ALIVE')
                        sj = mul(sig, sj)
                    cls[o_] = has_alive
                    ok_all &= has_alive
                per_sigma[(kind, 'out1 alive' if cls[0] else 'out1 none', 'out2 alive' if cls[1] else 'out2 none')] += 1
            print("conf %d prof %s lam %s %s: %d sigma (%.0fs)" % (ci, cf['prof'], cf['lam'], kind, n, time.time() - t0), flush=True)
    print("twisted producers (kind, output, status):")
    for k, v in sorted(tot.items()):
        print("   ", k, v)
    print("sigmas by whether each output's own-branch rows have an ALIVE non-class producer in the orbit:")
    for k, v in sorted(per_sigma.items()):
        print("   ", k, v)


if __name__ == "__main__":
    main()
