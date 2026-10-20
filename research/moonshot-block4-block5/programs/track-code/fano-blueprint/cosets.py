"""cosets.py -- ms-fano: the 16 cosets of RM(1,3) = extended Hamming [8,4,4] in F_2^8 as the frustration classes of
sign patterns on the 8 points of F_2^3 (pieces x1 + lam_v Z^v x2).  For each coset: leader weight, whether frustrated
(separable), smallest support of a frustrated sub-pattern, and the denominators of the exact separation (pattern.py)."""
import itertools
from pattern import separable, solve, check
pts = list(itertools.product((0, 1), repeat=3))
def aff(a, c):
    return tuple((c + sum(x * y for x, y in zip(a, v))) % 2 for v in pts)
code = {aff(a, c) for a in itertools.product((0, 1), repeat=3) for c in (0, 1)}
assert len(code) == 16
seen = set(); rows = []
for w in itertools.product((0, 1), repeat=8):
    if w in seen:
        continue
    coset = {tuple((x + y) % 2 for x, y in zip(w, c)) for c in code}
    seen |= coset
    leader = min(coset, key=sum)
    pat = [(v, -1 if b else 1) for v, b in zip(pts, leader)]
    ok, bad = separable(pat, 3)
    dens = None
    if ok:
        good, dens = check(pat, 3, solve(pat, 3)); assert good
    # smallest frustrated sub-pattern (support size) inside this leader's pattern
    best = None
    for r in range(2, 9):
        for S in itertools.combinations(range(8), r):
            sub = [pat[i] for i in S]
            if separable(sub, 3)[0]:
                best = r; break
        if best:
            break
    rows.append((sum(leader), leader, ok, best, dens))
rows.sort()
for wt, leader, ok, best, dens in rows:
    print(f"coset leader weight {wt}  signs {''.join('-' if b else '+' for b in leader)}  frustrated {ok}  "
          f"smallest frustrated sub-pattern {best}  denominators {dens}")
print(f"{len(rows)} cosets; covering radius {max(r[0] for r in rows)}")
