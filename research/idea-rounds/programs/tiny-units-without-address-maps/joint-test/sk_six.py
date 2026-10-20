#!/usr/bin/python3
"""joint/skeptic/sk_six.py -- agent "skeptic": do the six static inequalities alone exclude 3 arrays in 8 moves?
Rows taken from out_lp.txt (my model).  Also the two-line hand argument, checked by enumeration."""
import itertools
rows = [([2,4,4,4,4,2,0],7),([4,2,4,4,2,4,0],7),([4,4,2,2,4,4,0],7),([6,2,4,2,4,2,4],9),([2,6,4,2,2,4,4],9),([2,2,2,6,4,4,4],9)]
for W,k in ((3,8),(3,9),(2,5),(2,6),(1,2),(1,3),(5,13),(6,16)):
    tot=ok=0
    for c in itertools.combinations_with_replacement(range(7),k):
        n=[0]*7
        for j in c: n[j]+=1
        tot+=1
        if all(sum(a*b for a,b in zip(r,n))>=W*d for r,d in rows): ok+=1
    print("W=%d k=%d: multisets %d, satisfy all six inequalities: %d"%(W,k,tot,ok))
print("sum of the three pair rows:",[sum(r[0][j] for r in rows[:3]) for j in range(7)],"; sum of the three triple rows:",[sum(r[0][j] for r in rows[3:]) for j in range(7)])
print("all coefficients even:",all(a%2==0 for r,_ in rows for a in r))
