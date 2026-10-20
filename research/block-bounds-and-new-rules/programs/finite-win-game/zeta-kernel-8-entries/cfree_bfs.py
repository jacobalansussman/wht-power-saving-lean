# Cancellation-free in-place words for Z = T^(x)b, T = [[1,1],[0,1]] (row S = sum_{U >= S} x_U).
# A word is cancellation-free if no call ever cancels an entry; then the actual support of every register is
# its structural knowledge set K (union process), which only grows and must END at the target support.
# Register c starts with x_c, so c in K_c(end) = supp Z[pi(c)] = {U >= pi(c)}: pi(c) <= c for all c, which forces
# pi = id (induction on |c|).  So: BFS over tuples (K_c), K_c subset of U(c) = {U >= c}, calls i <- j allowed iff
# K_i u K_j subset of U(i).  Minimum number of calls to reach K_c = U(c) for all c.
import sys
from collections import deque

def run(b):
    N = 1 << b
    Ufull = [sum(1 << U for U in range(N) if (c & ~U) == 0) for c in range(N)]
    start = tuple(1 << c for c in range(N))
    goal = tuple(Ufull)
    dist = {start: 0}
    q = deque([start])
    while q:
        s = q.popleft()
        d = dist[s]
        if s == goal:
            return d, len(dist)
        for i in range(N):
            for j in range(N):
                if i == j: continue
                nk = s[i] | s[j]
                if nk == s[i] or (nk & ~Ufull[i]): continue
                t = s[:i] + (nk,) + s[i + 1:]
                if t not in dist:
                    dist[t] = d + 1
                    q.append(t)
    return None, len(dist)

if __name__ == "__main__":
    for b in (2, 3):
        d, n = run(b)
        print(f"b={b}: minimum cancellation-free in-place calls = {d}  (standard b*2^(b-1) = {b * 2 ** (b - 1)}; states explored {n})")
