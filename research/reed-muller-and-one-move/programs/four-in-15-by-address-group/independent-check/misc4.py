import sys
sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
GL = all_GL(); O4 = [g for g in GL if mul(g, tp(g)) == ID]
inv = {g: next(h for h in GL if mul(g, h) == ID) for g in O4}
def cls_of(s): return {frozenset(mul(mul(o, h), inv[o]) for h in s) for o in O4}
def gmul(a, b):
    r = 0
    for i in range(4):
        if (b >> i) & 1: r ^= a << i
    for d in (6, 5, 4):
        if (r >> d) & 1: r ^= 0b10011 << (d - 4)
    return r
def singer(k):
    ak = 1
    for _ in range(k): ak = gmul(ak, 2)
    cols = [gmul(ak, 1 << j) for j in range(4)]
    return frozenset(cyclic(tuple(tuple((cols[j] >> i) & 1 for j in range(4)) for i in range(4))))
reps = {15: ["4,9,7,12", "3,14,7,1", "3,15,8,14", "3,15,13,6", "3,4,10,5", "2,4,13,12", "3,9,14,15", "2,13,6,1", "3,4,14,11", "2,13,14,6"],
        5: ["3,8,13,15", "3,9,14,12", "2,10,14,15", "3,10,7,5", "4,10,11,6", "5,15,14,3", "3,6,10,13", "5,10,11,13", "2,6,10,3", "4,12,7,2"],
        3: ["6,8,15,10", "2,3,13,6", "6,13,10,3", "12,7,14,3", "4,8,5,10"]}
for k, n in [(1, 15), (3, 5), (5, 3)]:
    c = cls_of(singer(k))
    print('Singer power', k, '(order %d)' % n, 'lies in class of', [r for r in reps[n] if frozenset(cyclic(from_rows([int(v) for v in r.split(',')]))) in c])
