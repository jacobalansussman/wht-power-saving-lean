# Check 3: Alman-Li Lemma 1.1 (depth-2 exponent c  ->  depth-2d size d*2^{(1+c/d)k}) can only give Theta(N log N):
# min over integer d of d*2^{ck/d} / k  ->  c*e*ln2 for every fixed c > 0. So o(N log N) circuits (OpenAI; ours)
# cannot come from depth-2 circuits via this route unless the depth-2 exponent is 1 + o(1).
import math
for c in [0.5, 0.443, 0.1, 0.01]:
    row=[]
    for k in [20, 100, 1000, 10000]:
        best=math.exp(min(math.log(d)+(c*k/d)*math.log(2) for d in range(1,k+1)))/k
        row.append("k=%d: %.4f" % (k,best))
    print("c=%.3f  min_d d*2^(ck/d)/k:  %s   limit c*e*ln2 = %.4f" % (c, "  ".join(row), c*math.e*math.log(2)))
# Check 4: the same unordered semiring N[(Z>=1, x)] with 'Dirichlet-polynomial positivity on an interval' preorders:
# the inequality [3]+[1] >= 2[2] holds for n -> n^t iff 3^t + 1 >= 2*2^t; print where it holds on [0, 3].
ts=[i/20 for i in range(0,61)]
hold=[t for t in ts if 3**t+1 >= 2*2**t - 1e-12]
print("3^t+1 >= 2*2^t holds at t in", hold[:3], "...", "and fails on (0,1):", [t for t in ts if t not in hold][:5], "...")
