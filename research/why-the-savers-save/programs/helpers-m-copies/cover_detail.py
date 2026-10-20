"""print, for one group, the maximal giver chains (frames, givers, packets covered) - companion of chain_cover.py"""
import sys, json
exec(open('chain_cover.py').read().split("# best chains per final reader")[0].replace("gi = int(sys.argv[2])", "gi = int(sys.argv[2])"))
ports = w['ports']
for t in T:
    FL = [U for U in frames_list if contained(U, Ut[t])]
    succ = {U: [V for V in FL if len(V) > len(U) and contained(U, V)] for U in FL}
    res = []
    def dfs(U, m, path):
        m |= fr[U]; path = path + [U]
        res.append((m, path))
        for V in succ[U]: dfs(V, m, path)
    for U in FL: dfs(U, 0, [])
    res.sort(key=lambda x: -bin(x[0]).count('1'))
    seen = set()
    for m, path in res[:40]:
        if m in seen: continue
        seen.add(m)
        pc_ = Counter(pk_of[needed[k]] for k in range(len(needed)) if m >> k & 1)
        if bin(m).count('1') > 8:
            print('reader', t, 'cover', bin(m).count('1'), dict(pc_), 'path dims', [len(U) for U in path], 'givers', [sorted(frg[U]) for U in path])
