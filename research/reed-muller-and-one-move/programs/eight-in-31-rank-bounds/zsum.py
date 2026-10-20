"""zsum.py: summarise linkrec output: per subset, number of zero cells and zero-diagonal moves."""
import sys
for ln in open(sys.argv[1]):
    t = ln.split()
    if t[0] == 'SET': print(ln.strip()); continue
    if t[0] != 'Z': continue
    s = int(t[2]); Z = [int(h, 16) for h in t[3:]]
    nz = sum(bin(z).count('1') for z in Z)
    diag = [x + 1 for x in range(31) if Z[x] >> x & 1]
    print("S=%s |S|=%d zeros=%d diag=%d %s" % (format(s, '08b')[::-1], bin(s).count('1'), nz, len(diag), diag if len(diag) < 31 else 'all'))
