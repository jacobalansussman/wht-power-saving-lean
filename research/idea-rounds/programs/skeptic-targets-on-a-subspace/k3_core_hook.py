"""k3_core_hook.py (skeptic-2): Theorem R closes every no-scratch, non-interleaved 7-in-31 whose outside gates span
at most a 4-dimensional space of directions.

Set-up (e7-rm's reflection form, no scratch, 31 gates = each nonzero direction of F_2^5 once, 7 distinct invertible
targets T_l). O = outside gates, the other 31 - |O| gates are inside (sign flips of one array), all before and/or
after the outside ones. Then the outside part delivers C_l = T_l + I_l (I_l = sum of z z^T over array l's inside
gates), C_l lies in Sym(span O), and rank T_l = 5 <= rank C_l + n_l. Summing: sum_l rank C_l >= 35 - (31 - |O|)
= |O| + 4.  [argued in NOTE.md section 2; elementary]

Theorem R (L1) on the outside part, read as a program on span O: |O| * cap >= sum_l demand(C_l), and the hook demand
of a symmetric target depends only on its frame distance = its rank (P-double cosets).
Hook values: block 4, f = 9 (L9; cert-07 hookrel.py, v-o3-ratio local blocks, ms-fano mod p all give 68,320):
  psi = 27,456 / 48,576 / 68,320 / 86,592 at distance 1..4, cap 27,456.
Block 3, f = 5 (L7, Lean): 80 / 144 / 224, cap 80.
"""
import itertools
hooks = {4: ({0: 0, 1: 27456, 2: 48576, 3: 68320, 4: 86592}, 27456),
         3: ({0: 0, 1: 80, 2: 144, 3: 224}, 80)}
W = 7
for d, (psi, cap) in hooks.items():
    maxO = 2 ** d - 1
    for O in range(7, maxO + 1):
        need = O + 4                                   # sum of ranks the outside part must deliver
        best = None
        for ranks in itertools.combinations_with_replacement(range(d + 1), W):
            if sum(ranks) < need: continue
            dem = sum(psi[r] for r in ranks)
            if best is None or dem < best[0]: best = (dem, ranks)
        if best is None:
            print(f"span dim {d}, |O| = {O:2d}: rank {need} impossible with ranks <= {d} (max {W*d})"); continue
        margin = best[0] - O * cap
        print(f"span dim {d}, |O| = {O:2d}: need sum rank >= {need}; least demand {best[0]:>7,} {best[1]} "
              f"vs {O} x cap = {O*cap:>7,}  ->  {'EXCLUDED' if margin > 0 else 'open'} (margin {margin:,})")
# the cores the proposal names, as 7-array cores (4 arrays at rank 4, three arrays at 0)
print("four-alternating core / N17-type 4-in-15 core: sum of ranks =", 4 * 4, "< 19 required -> not a 7-array core")
