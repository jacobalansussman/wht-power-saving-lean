import sys
for f in sys.argv[1:]:
    L = open(f).read().split('\n'); vs = L[0]; p = L[1]; body = '\n'.join(L[2:]).strip().rstrip(',')
    open(f[:-3] + '.sing', 'w').write(f'ring r = {p}, ({vs}), dp;\noption(redSB);\nideal I = {body};\nint t = timer;\nideal G = std(I);\n'
        f'if (G[1] == 1) {{ print("GB=[1]"); }} else {{ print("CONSISTENT size " + string(size(G)) + " dim " + string(dim(G))); }}\nquit;\n')
