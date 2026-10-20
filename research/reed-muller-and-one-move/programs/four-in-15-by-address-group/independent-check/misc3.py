import sys
sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
GL = all_GL(); O4 = [g for g in GL if mul(g, tp(g)) == ID]
inv = {g: next(h for h in GL if mul(g, h) == ID) for g in O4}
def cls(rows):
    s = frozenset(cyclic(from_rows(rows)))
    return {frozenset(mul(mul(o, h), inv[o]) for h in s) for o in O4}
A = ["4,9,7,12", "3,14,7,1", "3,15,8,14", "3,15,13,6", "3,4,10,5", "2,4,13,12", "3,9,14,15", "2,13,6,1", "3,4,14,11", "2,13,14,6"]
B = ["6,4,13,10", "13,7,6,2", "3,11,10,14", "7,11,14,6", "9,5,10,4", "4,1,14,12", "11,13,12,14", "10,5,6,2", "9,13,6,12", "2,13,14,6"]
CA = [cls([int(v) for v in a.split(',')]) for a in A]
for i, b in enumerate(B):
    sb = frozenset(cyclic(from_rows([int(v) for v in b.split(',')])))
    print('clm P15%s (%s) lies in out_placements class %s' % ('abcdefghij'[i], b, ['abcdefghij'[j] for j in range(10) if sb in CA[j]]))
