import sys
sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
O3 = ["6,8,15,10", "2,3,13,6", "6,13,10,3", "12,7,14,3", "4,8,5,10"]
O5 = ["3,8,13,15", "3,9,14,12", "2,10,14,15", "3,10,7,5", "4,10,11,6", "5,15,14,3", "3,6,10,13", "5,10,11,13", "2,6,10,3", "4,12,7,2"]
O15 = ["4,9,7,12", "3,14,7,1", "3,15,8,14", "3,15,13,6", "3,4,10,5", "2,4,13,12", "3,9,14,15", "2,13,6,1", "3,4,14,11", "2,13,14,6"]
which = sys.stdin.readline().strip() if False else None
mode = open('sound_mode.txt').read().strip()
if mode == 'sound':
    npass = ntot = 0
    for rows in O3 + O5:
        H = cyclic(from_rows([int(v) for v in rows.split(',')]))
        for W, dl in [(1, [1, 2, 4, 8]), (2, [1, 2, 4, 8] * 2), (4, [1, 2, 4, 8] * 4), (4, [7, 11, 13, 14] * 4)]:
            g = Game(H, W); r = g.verdict(g.counts_from_dirs(dl)); ntot += 1; npass += r[0] == 'pass'
            print(rows, W, 'in', len(dl), dl[:4], r, flush=True)
    print('soundness: %d of %d known units pass' % (npass, ntot))
else:
    for rows in O15:
        H = cyclic(from_rows([int(v) for v in rows.split(',')]))
        for W, M in [(1, 3), (2, 7), (4, 15)]:
            g = Game(H, W); assert len(g.orbits) == 1
            print('P15', 'abcdefghij'[O15.index(rows)], rows, '%d in %d' % (W, M), g.verdict(g.counts_from_split((M,))), flush=True)
