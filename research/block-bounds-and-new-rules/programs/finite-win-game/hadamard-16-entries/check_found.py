# check_found.py -- deeper necessary checks for every zero-violation word the annealer printed
# (lines "g=.. restart .. best 0 : u-v u-v ..." in out_sa.txt):
#  (1) chk: every row set / column set of size <= 4 against all column / row sets (exact, generic mod p);
#  (2) the orthogonality lemma at EVERY downset of the word's heap (not only along its own order);
#  (3) the cell condition with future blocks of size 2 at every downset (both inputs know >= 8).
# Memory: downsets are stored as Python ints in a set; capped at 2e6 (reported if hit).
import re, subprocess, sys
src = sys.argv[1] if len(sys.argv) > 1 else 'out_sa.txt'
words = []
for line in open(src):
    m = re.match(r'g=(\d+) restart \d+ best 0 :(.*)', line.strip())
    if m:
        words.append([tuple(map(int, p.split('-'))) for p in m.group(2).split()])
def comps(calls, idx):
    lab = list(range(16))
    def f(a):
        while lab[a] != a: lab[a] = lab[lab[a]]; a = lab[a]
        return a
    for k in idx:
        u, v = calls[k]; a, b = f(u), f(v)
        if a != b: lab[a] = b
    d = {}
    for w in range(16): d.setdefault(f(w), set()).add(w)
    return list(d.values())
def downset_check(calls, cap=2_000_000):
    g = len(calls); last = {}; pred = []
    for k, (u, v) in enumerate(calls):
        pred.append([last[w] for w in (u, v) if w in last]); last[u] = k; last[v] = k
    seen = {0}; stack = [0]; bad = 0; single = 0; example = None
    while stack:
        D = stack.pop()
        inD = [k for k in range(g) if D >> k & 1]; out = [k for k in range(g) if not D >> k & 1]
        P = comps(calls, inD); Q = comps(calls, out)
        if not all(16 * len(A & B) == len(A) * len(B) for A in Q for B in P):
            bad += 1
            if example is None: example = (sorted(inD), [sorted(x) for x in P], [sorted(x) for x in Q])
        if len(P) * len(Q) == 16 and 2 <= len(P) <= 8: single += 1
        for k in out:
            if all(D >> j & 1 for j in pred[k]):
                E = D | (1 << k)
                if E not in seen:
                    seen.add(E); stack.append(E)
                    if len(seen) > cap: return len(seen), bad, single, example, True
    return len(seen), bad, single, example, False
for n, calls in enumerate(words):
    g = len(calls)
    print(f"--- word {n}: g={g}: " + ' '.join(f'{u}-{v}' for u, v in calls))
    inp = f"{g}\n" + '\n'.join(f'{u} {v}' for u, v in calls) + '\n'
    out = subprocess.run(['./chk'], input=inp, capture_output=True, text=True).stdout
    print('  chk:', out.strip().splitlines()[-1])
    nd, bad, single, ex, capped = downset_check(calls)
    print(f"  downsets {nd}{' (CAPPED)' if capped else ''}: orthogonality violated at {bad}; singleton-cell states {single}")
    if ex: print(f"  first violating downset: calls {ex[0]} past {ex[1]} future {ex[2]}")
    sys.stdout.flush()
