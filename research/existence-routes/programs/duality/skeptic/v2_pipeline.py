# v-ex2-duality session 2: literal simulation of the Pipelining Lemma's schedule (ex2-duality PROOF.md section 2).
# Q makes L sequential calls (call j needs the results of calls 1..j-1 of the same copy); N copies are staggered;
# in round t copy i issues call t-i+1.  P serves a set of pending calls only if it is one call of EACH number 1..L
# (one copy of y/N, given simultaneously).  Other rounds ("boundary") are served directly (plain + A5).
# Checks: dependency order respected; every call served once; steady rounds = N-L+1; boundary calls = L(L-1).
import random
def run(L, N):
    done = {}            # (copy, call) -> round served
    steady, boundary_calls = 0, 0
    for t in range(1, N + L):
        pend = [(i, t - i + 1) for i in range(1, N + 1) if 1 <= t - i + 1 <= L]
        for (i, j) in pend:   # dependency: calls 1..j-1 of copy i served in earlier rounds
            assert all((i, jj) in done and done[(i, jj)] < t for jj in range(1, j)), (i, j, t)
        nums = sorted(j for _, j in pend)
        if nums == list(range(1, L + 1)): steady += 1
        else: boundary_calls += len(pend)
        for p in pend:
            assert p not in done; done[p] = t
    assert len(done) == N * L
    return steady, boundary_calls
ok = True
for _ in range(300):
    L, N = random.randint(1, 12), random.randint(1, 60)
    s, b = run(L, N)
    exp_s = max(0, N - L + 1)
    exp_b = N * L - L * exp_s
    if s != exp_s or b != exp_b: ok = False; print("MISMATCH", L, N, s, b)
    if N >= L and b != L * (L - 1): ok = False; print("boundary count", L, N, b)
print("pipeline schedule: dependency order OK in all runs; steady rounds = N-L+1 and boundary calls = L(L-1) when N>=L:", ok)
print("so the boundary overhead is L(L-1) calls of weight 1/N each (O(L^2/N), not O(L/N) as PROOF.md writes; still -> 0 as N grows)")
