# Independent re-pricing of C588MHH's block list (from out/C588MHH.json) with the Lean-route price pipeline's price() (exact root + plan_block at s=40)
import sys, json
sys.path.insert(0, __import__('os').path.dirname(__file__))
from route_price import price
d=json.load(open(sys.argv[1]))
H={int(k):v for k,v in d['hist'].items()}
W=d['W']; m=73
D=m*W-sum(r*n for r,n in H.items())
print('W',W,'D',D,'mass',sum(r*n for r,n in H.items()), 'count', sum(H.values()))
print(price(dict(hist=H,W=W,D=D)))
