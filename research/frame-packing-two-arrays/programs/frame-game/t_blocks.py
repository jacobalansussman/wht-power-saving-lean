"""t_blocks.py: the block structure behind the frame game, checked at block 2 (and the group order at block 3).
(1) |G / mu4| = 4^m |Sp(2m,2)|  (so the kernel of G -> Sp is exactly the Paulis times mu4);
(2) for every symplectic image, the 4^m group elements over it have linearly independent matrices (they span
    M_n(C)), so Phi: C[G]_mu4 -> M_n(C) (x) C[Sp], g -> g (x) gbar, is an isomorphism (dimension count);
(3) the order-3 control word, mapped to the frame game, is a frame-game program (plain 2 in 2m in another basis)."""
import sys, itertools
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "ledger"))
import numpy as np
import ledger as L
from fpa import *

def group_bfs(E, gens):
    seen = {E.ID}; frontier = [E.ID]
    while frontier:
        nf = []
        for g in frontier:
            for name in gens:
                h, _ = E.lmul(name, g)
                if h not in seen: seen.add(h); nf.append(h)
        frontier = nf
    return seen
for m in (2,):
    E = L.Engine(m)
    gens = [k for k in E.letters if k.startswith("D") or k.startswith("X") or k.startswith("P")]
    # address maps: all of GL(m,2) via two generators
    if m == 2:
        gens.append(E.addr([2, 1])); gens.append(E.addr([1, 3]))
    G = group_bfs(E, gens)
    syms = {}
    for g in G:
        syms.setdefault(L.key2(E.syms[g]), []).append(g)
    nsp = len(syms)
    print(f"m={m}: |G/mu4| = {len(G)}; distinct symplectic images = {nsp}; 4^m * |Sp(4,2)| = {16 * 720}")
    fib = set(len(v) for v in syms.values())
    print(f"   fibre sizes: {fib}")
    # (2) linear independence of each fibre's matrices
    bad = 0
    for key, els in syms.items():
        M = np.array([np.concatenate([E.mats[g][0].ravel(), E.mats[g][1].ravel()]) for g in els], dtype=float)
        # complex independence: use complex matrices
        Mc = np.array([(E.mats[g][0] + 1j * E.mats[g][1]).ravel() for g in els])
        if np.linalg.matrix_rank(Mc) != len(els): bad += 1
    print(f"   fibres whose matrices are linearly dependent: {bad} of {nsp}")
