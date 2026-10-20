import sys
sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
GL = all_GL(); O4 = [g for g in GL if mul(g, tp(g)) == ID]
inv = {g: next(h for h in GL if mul(g, h) == ID) for g in O4}
def cls(rows):
    s = frozenset(cyclic(from_rows(rows)))
    return {frozenset(mul(mul(o, h), inv[o]) for h in s) for o in O4}
e3 = cls([12, 5, 11, 6]); rep = frozenset(cyclic(from_rows([6, 13, 10, 3])))
print('alpha^5 class size', len(e3), 'contains rep [6,13,10,3]:', rep in e3)
# symmetric generators per order-15 class
O15 = ["4,9,7,12", "3,14,7,1", "3,15,8,14", "3,15,13,6", "3,4,10,5", "2,4,13,12", "3,9,14,15", "2,13,6,1", "3,4,14,11", "2,13,14,6"]
for i, r in enumerate(O15):
    c = cls([int(v) for v in r.split(',')])
    sym = any(any(h == tp(h) and len(cyclic(h)) == 15 for h in s) for s in c)
    print('P15' + 'abcdefghij'[i], r, 'class size', len(c), 'has symmetric generator:', sym)
for r in ["6,8,15,10", "2,3,13,6", "6,13,10,3", "12,7,14,3", "4,8,5,10"]:
    c = cls([int(v) for v in r.split(',')])
    print('order3', r, len(c), 'symmetric generator:', any(any(h == tp(h) and h != ID for h in s) for s in c))
for r in ["3,8,13,15", "3,9,14,12", "2,10,14,15", "3,10,7,5", "4,10,11,6", "5,15,14,3", "3,6,10,13", "5,10,11,13", "2,6,10,3", "4,12,7,2"]:
    c = cls([int(v) for v in r.split(',')])
    print('order5', r, len(c), 'symmetric generator:', any(any(h == tp(h) and h != ID for h in s) for s in c))
# 1 in 3 at orders 3 and 5 (one-orbit-free splits: all splits of 3)
for r in ["6,8,15,10", "2,3,13,6", "6,13,10,3", "12,7,14,3", "4,8,5,10", "3,8,13,15", "3,9,14,12", "2,10,14,15", "3,10,7,5", "4,10,11,6", "5,15,14,3", "3,6,10,13", "5,10,11,13", "2,6,10,3", "4,12,7,2"]:
    g = Game(cyclic(from_rows([int(v) for v in r.split(',')])), 1)
    res = [g.verdict(g.counts_from_split(sp))[0] for sp in splits(3, len(g.orbits))]
    print('1 in 3', r, 'passing splits:', sum(x == 'pass' for x in res), 'of', len(res))
