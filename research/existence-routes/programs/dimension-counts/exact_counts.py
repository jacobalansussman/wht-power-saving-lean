# Exact-level (group-ring) counts for the normal form of pf-algebraic Theorem 8: unknowns a_{j,z}, f_{t,z} in C[B];
# contents in C[G]^W.  All quantities in units of |B|.  Group orders from Lemma 6 (|G| = 4 * 4^m * |Sp(2m,2)|,
# |B| = 4 * 4^m * |P|), Lagrangians at distance k from L0: n_k = 2^(k(k+1)/2) [m k]_2, L = sum_k n_k = [G:B].
from fractions import Fraction as Fr
from math import prod
def gauss_binom(m, k, q=2):
    return prod(q**(m - i) - 1 for i in range(k)) // prod(q**(i + 1) - 1 for i in range(k))
def nk(m, k): return 2**(k*(k+1)//2) * gauss_binom(m, k)
def sp_order(m): return 2**(m*m) * prod(4**i - 1 for i in range(1, m+1))
def gl_order(m): return prod(2**m - 2**i for i in range(m))
def counts(m, W, steps, targets):
    n = len(steps); L = sum(nk(m, k) for k in range(m+1))
    slots = sum(W + j for j in range(n)) + W * (W + n)
    E = W * W * L
    # orbit of the target tuple under post-mixing by M_W(C[B]) and pre-mixing by M_W(C[B]):
    # entry (t,s): dim(C[B] T_s + T_t C[B]) = 2|B| - |B T_s cap T_t B|; = 2 - 1/n_k if T_s, T_t in the same double coset
    o = sum(2 - (Fr(1, nk(m, targets[t])) if targets[t] == targets[s] else 0) for t in range(W) for s in range(W))
    gauge = sum(Fr(1, nk(m, r)) for r in steps)          # step stabilisers B cap W_r^-1 B W_r (a lower bound on the fibre)
    return dict(L=L, slots=slots, P=slots, E=E, naive=slots - E, orbit=float(o), gauge=float(gauge),
                Pmax_minus_gauge=float(slots - gauge), orbit_test_passes=float(slots - gauge) >= float(o),
                dominance_possible=slots - float(gauge) >= E)
shapes = [
 # label, m, W, steps, targets, status
 ('block 2 plain 1 in 2', 2, 1, [1,1], [2], 'exists'),
 ('block 2 plain 2 in 4', 2, 2, [1]*4, [2,2], 'exists'),
 ('block 2 1 in 1', 2, 1, [1], [2], 'empty (proof)'),
 ('block 2 2 in 3', 2, 2, [1]*3, [2,2], 'empty (proof)'),
 ('block 2 3 in 5 (XW target)', 2, 3, [1]*5, [2,2,2], 'engine: empty; XW: open'),
 ('block 3 plain 1 in 3', 3, 1, [1]*3, [3], 'exists'),
 ('block 3 plain 2 in 6', 3, 2, [1]*6, [3,3], 'exists'),
 ('block 3 2 in 5', 3, 2, [1]*5, [3,3], 'empty (proof); XW: open'),
 ('block 3 3 in 8', 3, 3, [1]*8, [3,3,3], 'empty (14W/5); XW: open'),
 ('block 3 5 in 14', 3, 5, [1]*14, [3]*5, 'open'),
 ('block 3 mixed 1+4 in 10', 3, 5, [1]*10, [3,2,2,2,2], 'open'),
 ('block 4 plain 1 in 4', 4, 1, [1]*4, [4], 'exists'),
 ('block 4 plain 2 in 8', 4, 2, [1]*8, [4,4], 'exists'),
 ('block 4 2 in 7', 4, 2, [1]*7, [4,4], 'open'),
 ('block 4 2 in 7 grouped (2,2,2,1)', 4, 2, [2,2,2,1], [4,4], 'open'),
 ('block 4 3 in 10', 4, 3, [1]*10, [4,4,4], 'open'),
]
if __name__ == '__main__':
    print(f"{'shape':36s} {'status':26s} {'L':>5s} {'P/|B|':>6s} {'E/|B|':>7s} {'naive':>7s} {'orbit':>7s} {'P-gauge':>8s} dominance-possible")
    for lab, m, W, st, tg, status in shapes:
        c = counts(m, W, st, tg)
        print(f"{lab:36s} {status:26s} {c['L']:5d} {c['P']:6d} {c['E']:7d} {c['naive']:7d} {c['orbit']:7.3f} {c['Pmax_minus_gauge']:8.3f} {c['dominance_possible']}")
