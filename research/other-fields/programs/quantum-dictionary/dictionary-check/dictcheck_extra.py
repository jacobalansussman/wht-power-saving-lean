#!/usr/bin/env python3
"""x-dict extra checks (13 Oct 2026): (H) OpenAI family-130 frames Phi_U = J diag(i^{wt(P_U x)}) J have label L_U;
(I) every Lagrangian = (X-projection U, symmetric form beta on U); subspace frames = beta = dot|U, CSS = beta = 0.
Run: python3 -I dictcheck_extra.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dictcheck import *


def section_H():
    say("== H. OpenAI family-130 frames (network.tex eq. net:frame): L(Phi_U) = L_U for nondegenerate U")
    for m in (1, 2, 3):
        n = 1 << m
        W = tensor_power(littleH(), m)
        cnt = 0
        for U in all_subspaces(m):
            Up = perp_dot(U, m)
            if len(U & Up) != 1:
                continue                     # degenerate: Phi_U undefined in the paper
            proj = {}
            for u in U:
                for v in Up:
                    proj[u ^ v] = u
            D = [[IPOW[popc(proj[x]) % 4] if x == y else Z0 for y in range(n)] for x in range(n)]
            Phi = mmul(mmul(W, D), W)        # 2^m Phi_U
            assert label_of(Phi, m) == L_sub(U, m)
            cnt += 1
        say("  m=%d: %d nondegenerate U, label of Phi_U = L_U = {(u|u+v)} for all: True" % (m, cnt))


def section_I():
    say("== I. Lagrangian = (U = projection to the X part, symmetric form beta on U); counts sum_k [m,k]_2 2^(k(k+1)/2)")
    for m in (1, 2, 3, 4):
        Ls = all_lagrangians_direct(m)
        M = (1 << m) - 1
        by = {}
        for L in Ls:
            U = frozenset(v & M for v in L)
            by.setdefault(U, []).append(L)
        ok_count = all(len(by[U]) == 2 ** (dim2(U) * (dim2(U) + 1) // 2) for U in by)
        # subspace frame over U is the one with beta = dot product: contains (u|u) for every u in U
        sub_ok = all(sum(1 for L in by[U] if all((u | (u << m)) in L for u in U)) == 1 and
                     [L for L in by[U] if all((u | (u << m)) in L for u in U)][0] == L_sub(U, m) for U in by)
        css_ok = all(sum(1 for L in by[U] if all(u in L for u in U)) == 1 for U in by)
        say("  m=%d: %d X-projections (all subspaces: %s); each carries 2^(k(k+1)/2) Lagrangians: %s; exactly one with"
            " beta = dot (the subspace frame L_U): %s; exactly one with beta = 0 (CSS): %s; non-subspace frames: %d"
            % (m, len(by), len(by) == len(all_subspaces(m)), ok_count, sub_ok, css_ok, len(Ls) - len(by)))
        assert ok_count and sub_ok and css_ok


if __name__ == "__main__":
    section_H()
    section_I()
    say("EXTRA CHECKS PASSED")
