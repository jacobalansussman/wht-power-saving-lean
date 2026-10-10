#!/usr/bin/env python3
"""groups.py -- generator sets (label maps) of the groups tried, all inside the design group O+(8,2)."""
import sys
sys.dont_write_bytecode = True
from lib import *
def f8_perms():
    """PSL(2,8) on the projective line over F_8 = F_2[a]/(a^3 + a + 1): points 0..7 = field elements, 8 = infinity"""
    def mul(x, y):
        r = 0
        for i in range(3):
            if y >> i & 1: r ^= x << i
        for i in (4, 3):
            if r >> i & 1: r ^= 0b1011 << (i - 3)
        return r
    inv = {x: next(y for y in range(1, 8) if mul(x, y) == 1) for x in range(1, 8)}
    tr = [x ^ 1 for x in range(8)] + [8]                 # x -> x + 1
    sc = [mul(2, x) for x in range(8)] + [8]             # x -> a x
    iv = [8] + [inv[x] for x in range(1, 8)] + [0]       # x -> 1/x
    fr = [mul(x, x) for x in range(8)] + [8]             # Frobenius
    return tr, sc, iv, fr
def group(name):
    tr, sc, iv, fr = f8_perms()
    if name == "trivial": return []
    if name == "C9": return [perm_map([1, 2, 3, 4, 5, 6, 7, 8, 0])]
    if name == "C3xC3": return [perm_map([1, 2, 0, 4, 5, 3, 7, 8, 6]), perm_map([3, 4, 5, 6, 7, 8, 0, 1, 2])]
    if name == "omega": return [omega()]
    if name == "C3": return [perm_map([1, 2, 0, 4, 5, 3, 7, 8, 6])]
    if name == "C3b": return [perm_map([3, 4, 5, 6, 7, 8, 0, 1, 2])]
    if name == "C3cubed": return [perm_map([1, 2, 0, 3, 4, 5, 6, 7, 8]), perm_map([0, 1, 2, 4, 5, 3, 6, 7, 8]), perm_map([0, 1, 2, 3, 4, 5, 7, 8, 6])]
    if name == "C3xC3om": return [perm_map([1, 2, 0, 4, 5, 3, 7, 8, 6]), perm_map([3, 4, 5, 6, 7, 8, 0, 1, 2]), omega()]
    if name == "C3om": return [perm_map([3, 4, 5, 6, 7, 8, 0, 1, 2]), omega()]
    if name == "AGL23": return [perm_map([1, 2, 0, 4, 5, 3, 7, 8, 6]), perm_map([3, 4, 5, 6, 7, 8, 0, 1, 2]), perm_map([0, 2, 1, 6, 8, 7, 3, 5, 4])]
    if name == "C3xC3_2": return [perm_map([1, 2, 0, 4, 5, 3, 7, 8, 6]), perm_map([3, 4, 5, 6, 7, 8, 0, 1, 2]), perm_map([0, 2, 1, 6, 8, 7, 3, 5, 4])][:3]
    if name == "C2": return [perm_map([1, 0, 3, 2, 5, 4, 7, 6, 8])]
    if name == "V4": return [perm_map([1, 0, 3, 2, 5, 4, 7, 6, 8]), perm_map([2, 3, 0, 1, 6, 7, 4, 5, 8])]
    if name == "C5": return [perm_map([1, 2, 3, 4, 0, 5, 6, 7, 8])]
    if name == "C4": return [perm_map([1, 2, 3, 0, 5, 6, 7, 4, 8])]
    if name == "C6": return [perm_map([1, 2, 0, 4, 5, 3, 7, 6, 8])]
    if name.startswith("orth"):
        # reflections in k mutually orthogonal roots (ports that pairwise do not know each other), found greedily
        k = int(name[4:]); ch = []
        for i in range(V):
            if all(i != j and not e8.knows(i, j) for j in ch):
                ch.append(i)
            if len(ch) == k: break
        return [transvect(P[i]) for i in ch]
    if name == "E16": return [perm_map([1, 0, 3, 2, 4, 5, 6, 7, 8]), perm_map([2, 3, 0, 1, 4, 5, 6, 7, 8]), perm_map([0, 1, 2, 3, 5, 4, 7, 6, 8]), perm_map([0, 1, 2, 3, 6, 7, 4, 5, 8])]
    if name == "D8": return [perm_map([1, 2, 3, 0, 4, 5, 6, 7, 8]), perm_map([1, 0, 3, 2, 4, 5, 6, 7, 8])]
    if name == "F8orth":
        tr, sc, iv, fr = f8_perms()
        return [perm_map(tr), perm_map([x ^ 2 for x in range(8)] + [8]), perm_map([x ^ 4 for x in range(8)] + [8]), transvect(0b011111111 ^ 0b100000000 ^ 0b111111111 ^ 0b011111111 | 0b111)]
    if name == "V4b": return [perm_map([1, 0, 3, 2, 4, 5, 6, 7, 8]), perm_map([2, 3, 0, 1, 4, 5, 6, 7, 8])]
    if name == "E8b": return [perm_map([1, 0, 2, 3, 4, 5, 6, 7, 8]), perm_map([0, 1, 3, 2, 4, 5, 6, 7, 8]), perm_map([0, 1, 2, 3, 5, 4, 6, 7, 8])]
    if name == "E16b": return [perm_map([1, 0, 2, 3, 4, 5, 6, 7, 8]), perm_map([0, 1, 3, 2, 4, 5, 6, 7, 8]), perm_map([0, 1, 2, 3, 5, 4, 6, 7, 8]), perm_map([0, 1, 2, 3, 4, 5, 7, 6, 8])]
    if name == "E32": return group("E16") + [perm_map([4, 5, 6, 7, 0, 1, 2, 3, 8])]
    if name == "Syl128": return [perm_map([1, 0, 2, 3, 4, 5, 6, 7, 8]), perm_map([2, 3, 0, 1, 4, 5, 6, 7, 8]), perm_map([4, 5, 6, 7, 0, 1, 2, 3, 8])]
    if name == "D8xD8": return [perm_map([1, 2, 3, 0, 4, 5, 6, 7, 8]), perm_map([2, 1, 0, 3, 4, 5, 6, 7, 8]), perm_map([0, 1, 2, 3, 5, 6, 7, 4, 8]), perm_map([0, 1, 2, 3, 6, 5, 4, 7, 8])]
    if name == "F8u": return group("F8plus") + [perm_map([x ^ ((x >> 1) & 1) for x in range(8)] + [8])]
    if name == "F8uu": return group("F8u") + [perm_map([x ^ (((x >> 2) & 1) << 1) for x in range(8)] + [8])]
    if name == "E16orth": return group("E16") + group("orth3")[:1]
    if name == "E16x": return [perm_map([1, 0, 3, 2, 4, 5, 6, 7, 8]), perm_map([2, 3, 0, 1, 4, 5, 6, 7, 8]), perm_map([0, 1, 2, 3, 5, 4, 7, 6, 8]), perm_map([0, 1, 2, 3, 6, 7, 4, 5, 8]), transvect(0b011111111 & ~0b1 | 0b100000000) if False else transvect(P[0])]
    if name == "E8c": return [perm_map([1, 0, 3, 2, 4, 5, 6, 7, 8]), perm_map([2, 3, 0, 1, 4, 5, 6, 7, 8]), perm_map([0, 1, 2, 3, 5, 4, 7, 6, 8])]
    if name == "E8d": return [perm_map([1, 0, 3, 2, 5, 4, 7, 6, 8]), perm_map([2, 3, 0, 1, 6, 7, 4, 5, 8]), perm_map([0, 1, 2, 3, 5, 4, 7, 6, 8])]
    if name == "SylR1": return group("Syl128") + [transvect(0b100000011)]
    if name == "SylR2": return group("Syl128") + [transvect(0b000000111)]
    if name == "SylR3": return group("Syl128") + [compose(transvect(0b000001111 ^ 0b1000), transvect(0b011110000 ^ 0b10000000))]
    if name == "F8uuu": return group("F8uu") + [perm_map([x ^ ((x >> 2) & 1) for x in range(8)] + [8])]
    if name == "refl": return [transvect(0b111)]
    if name == "refl2": return [transvect(0b111), transvect(0b111000)]
    if name == "C3xC3refl": return [perm_map([1, 2, 0, 4, 5, 3, 7, 8, 6]), perm_map([3, 4, 5, 6, 7, 8, 0, 1, 2]), compose(transvect(0b111), transvect(0b111000), transvect(0b111000000))]
    if name == "PSL28": return [perm_map(tr), perm_map(sc), perm_map(iv)]
    if name == "PGammaL28": return [perm_map(tr), perm_map(sc), perm_map(iv), perm_map(fr)]
    if name == "S8": return [perm_map([1, 0, 2, 3, 4, 5, 6, 7, 8]), perm_map([1, 2, 3, 4, 5, 6, 7, 0, 8])]
    if name == "S9": return [perm_map([1, 0, 2, 3, 4, 5, 6, 7, 8]), perm_map([1, 2, 3, 4, 5, 6, 7, 8, 0])]
    if name == "AGL18": return [perm_map(tr), perm_map(sc)]                  # order 56, fixes infinity
    if name == "F8plus": return [perm_map(tr), perm_map([x ^ 2 for x in range(8)] + [8]), perm_map([x ^ 4 for x in range(8)] + [8])]   # 2^3
    if name == "C7": return [perm_map(sc)]
    if name == "full": return [perm_map([1, 0, 2, 3, 4, 5, 6, 7, 8]), perm_map([1, 2, 3, 4, 5, 6, 7, 8, 0]), transvect(0b111)]
    if name in ("spread", "spreadN"):
        # complex reflections of the 40 lines of omega: products of the two reflections of a line; they commute with omega
        om = omega(); p = as_perm(om)
        gs = [om]
        for a in (0, 1, 5, 17, 40, 90):
            gs.append(compose(transvect(P[a]), transvect(P[p[a]])))
        return gs
    raise KeyError(name)
def order(gens, cap=200000):
    """order of the permutation group on the 120 ports (closure by BFS on permutations), None above cap"""
    perms = [tuple(as_perm(g)) for g in gens]
    ident = tuple(range(V)); seen = {ident}; todo = [ident]
    while todo:
        x = todo.pop()
        for p in perms:
            y = tuple(p[i] for i in x)
            if y not in seen:
                seen.add(y); todo.append(y)
                if len(seen) > cap: return None
    return len(seen)
if __name__ == "__main__":
    for nm in ("C9", "C3xC3", "omega", "C7", "F8plus", "AGL18", "PSL28", "PGammaL28", "spread"):
        gs = group(nm)
        ok = all(as_perm(g) is not None for g in gs)
        print(nm, "permute ports:", ok, " order:", order(gs, 80000))
