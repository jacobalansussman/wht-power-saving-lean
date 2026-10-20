#!/usr/bin/env python3
"""vfieldvec.py (v-rm-singer): vectorised second code path for S3/S4/S5 over ALL forms at m = 3..7.
Uses F_2-linearity in the form coefficients: the move vector K(Q) (bits k: Q(a^-k)) and the target vector J(Q)
(bits k: <Q, s^-k s^-kT>, from the matrix) are XORs of the basis-monomial vectors.  Stabiliser of Q in Z_n:
rotation of K by n/p for every prime p | n.  Also checks the rotation law a_{Q o s} = a * a_Q on random forms.
"""
import sys, itertools, random
import numpy as np
sys.path.insert(0, '.')
from vfield import GF, setup, channel_data, aQ, Qval


def popcount64(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return ((x * np.uint64(0x0101010101010101)) >> np.uint64(56)).astype(np.int64)


def main():
    ms = [int(x) for x in sys.argv[1].split(',')]
    for m in ms:
        for coords in ('sd', 'poly'):
            F, basis, C, sig, Sk, pairs = setup(m, coords)
            n = F.n; N = len(pairs)
            assert n <= 63
            bK = []; bJ = []
            for j in range(N):
                q = tuple(int(i == j) for i in range(N))
                K, J = channel_data(F, C, Sk, pairs, q)
                bK.append(sum(b << k for k, b in enumerate(K))); bJ.append(sum(b << k for k, b in enumerate(J)))
            # all 2^N forms by XOR combination
            idx = np.arange(1 << N, dtype=np.int64)
            Kv = np.zeros(1 << N, dtype=np.uint64); Jv = np.zeros(1 << N, dtype=np.uint64)
            for j in range(N):
                sel = ((idx >> j) & 1).astype(bool)
                Kv[sel] ^= np.uint64(bK[j]); Jv[sel] ^= np.uint64(bJ[j])
            wK = popcount64(Kv); wJ = popcount64(Jv)
            nz = idx != 0
            assert (wK[nz] > 0).all()
            odd = (wJ > 0)
            mask = np.uint64((1 << n) - 1)
            def rot(v, r):
                r = np.uint64(r)
                return ((v << r) | (v >> (np.uint64(n) - r))) & mask
            primes = [p for p in range(2, n + 1) if n % p == 0 and all(p % d for d in range(2, p))]
            stab_nontriv = np.zeros(1 << N, dtype=bool)
            for p in primes:
                stab_nontriv |= (rot(Kv, n // p) == Kv)
            stab_nontriv &= nz
            ratio = np.where(odd, wJ / np.maximum(wK, 1), 0.0)
            imax = int(np.argmax(ratio))
            line = (f"m={m} {coords:4s} n={n} N={N}: nonzero forms {nz.sum()}, odd {odd.sum()}"
                    f" | odd with nontrivial stabiliser: {(odd & stab_nontriv).sum()}"
                    f" | even with nontrivial stabiliser: {((~odd) & stab_nontriv).sum()}"
                    f" | min |K| over nonzero {wK[nz].min()}, over odd {wK[odd].min()} (2^(m-2)={2 ** (m - 2)})"
                    f" | odd |J| values {sorted(set(wJ[odd].tolist()))}"
                    f" | deletion max |J|/|K| = {ratio[imax]:.6f} at |K|={wK[imax]}, |J|={wJ[imax]}")
            if coords == 'sd':
                # orbit count of odd forms: all free -> odd/n
                line += f" | odd orbits = odd/n = {odd.sum() / n} (2^(N-m) = {2 ** (N - m)})"
            print(line, flush=True)
            if coords == 'sd':
                rng = random.Random(7 + m); bad = 0
                for _ in range(200):
                    q = tuple(rng.randint(0, 1) for _ in range(N))
                    # Q o s as a function: y -> Q(a y); recover a_{Q o s} directly
                    a1 = 0
                    for y in range(1, F.q):
                        if Qval(q, pairs, C.c(F.mul(F.pw(1), y))): a1 ^= F.inv(y)
                    if a1 != F.mul(F.pw(1), aQ(F, C, pairs, q)): bad += 1
                print(f"   rotation law a_(Q o s) = a a_Q: {bad} violations in 200 random forms", flush=True)


if __name__ == '__main__':
    main()
