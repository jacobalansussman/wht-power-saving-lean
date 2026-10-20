import vsurv as S
v = {1: -1, 2: 4, 3: -5, 4: 2}
res = S.tests(v, 4, Mextra=4)
bad = sorted([r for r in res if r[1] < 0], key=lambda r: r[1])
print(len(res), "tests;", len(bad), "fail; first:", [(n, int(x)) for n, x in bad[:12]])
v5 = {2: -1, 3: 4, 4: -5, 5: 2}
res = S.tests(v5, 5, Mextra=3); bad = [r for r in res if r[1] < 0]
print("[2]+5[4] <= 4[3]+2[5]:", len(res), "tests;", len(bad), "fail;", [(n, int(x)) for n, x in bad[:8]])
