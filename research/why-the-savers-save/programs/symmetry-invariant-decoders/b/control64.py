"""control64.py (t-km): control for km7b's exact tests: the hub's natural decoder on P112 (8 coordinate totals,
frames = spans of the ports non-orthogonal to each removed octad line, cost 64) must PASS stages 0-2.
Also checks whether that 8-frame family is invariant under g7 / g3."""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
sys.argv = ['km7b.py', '7', '8', '0', '99']
src = open('km7b.py').read().split("t0 = time.time()")[0]
exec(src)
import itertools
exec("def full_tests(frames):" + open('km7b.py').read().split("def full_tests(frames):")[1].split("stats = Counter()")[0])
frames = []
for o in list(O0):
    reads = [i for i, T in enumerate(P) if GRAM[o, T] != 0]
    b = f2basis([LAB[i] for i in reads]); sp = {0}
    for x in b: sp |= {y ^ x for y in sp}
    frames.append(np.array([LAB[i] in sp for i in range(NP)]))
print('natural frames: dims', [f2dim([LAB[i] for i in np.nonzero(f)[0]]) for f in frames], 'ports', [int(f.sum()) for f in frames])
print('full tests verdict on the natural decoder:', full_tests(frames))
for q in (7, 3):
    pm = np.load('gens.npz')['perm%d' % q]
    fs = {f.tobytes() for f in frames}
    print('natural family invariant under g%d:' % q, all(f[np.argsort(pm)].tobytes() in fs for f in frames))
