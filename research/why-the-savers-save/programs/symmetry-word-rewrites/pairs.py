import sys; sys.path.insert(0,'.')
import fw, search, reorder
from fractions import Fraction as Fr
from collections import Counter
W = search.Word(fw.load('K8-FULL-WORD.json.gz'))
O = reorder.OW(W)
seq = W.seq
st = Counter(); ex = []
N = len(seq)
for i in range(N):
    e = seq[i][2]
    if e[0] != 'add': continue
    for j in range(i+1, N):
        h = seq[j][2]
        if O.sec[j] != O.sec[i]: break
        if h[0]=='add' and h[1]==e[1] and h[2]==e[2]:
            # can they meet? all between commute with e (move e down) or with h (move h up)
            ok = all(O.commute(i, k) for k in range(i+1, j)) or all(O.commute(j, k) for k in range(i+1, j))
            q = Fr(e[3],e[4]) + Fr(h[3],h[4])
            st[(O.sec[i], 'adjacentable' if ok else 'blocked', 'cancel' if q==0 else 'merge')] += 1
            if ok and len(ex)<10: ex.append((i,j,e,h))
            break
        if not O.commute(i, j) and (h[1] in (e[1],e[2]) or (h[0]=='add' and h[2] in (e[1],))):
            # stop at first conflict involving same registers
            pass
print(st); print(ex[:10])
