import sys, json, math
sys.path.insert(0, __import__('os').path.dirname(__file__))
from route_price import unit, price, widths_from, IDLE
from exact_moment import certify
from collections import Counter
m=73
def L(h): return sum(n*r*math.log(m/r) for r,n in h.items())
def root(h,W): return float(certify({int(k):v for k,v in h.items()},m,W)['root_lower_decimal'][:18])
nat588={1:173404,2:40632,3:23928,4:14024,5:4038,6:4032,8:576,30:1,48:1857}
print('native #588 hosted (hub, W=7168): root', root(nat588,7168), ' L=%.0f'%L(nat588), ' 2560/L=%.6e'%(2560/L(nat588)))
loc={1:2484,2:492,3:280,4:232,5:56,6:56,8:8}; wid={4:48,5:8,7:32,8:140,9:276}
u0=unit(loc,widths_from(wid),copies=False)
print('banks, NO source copies: W',u0['W'],'D',u0['D'],'root',root(u0['hist'],u0['W']), 'L=%.0f'%L(u0['hist']))
u1=unit(loc,widths_from(wid),copies=True)
print('banks + source copies (Shea): W',u1['W'],'D',u1['D'],'root',root(u1['hist'],u1['W']),'L=%.0f'%L(u1['hist']))
# L by component for Shea
comp={}
comp['72 x local word']=L({r:72*n for r,n in loc.items()})
comp['72 x source copies (rank1+rank8 per source)']=72*112*(1*math.log(73)+8*math.log(73/8))
comp['idle charge of the 112 source banks {8:48,16:32,32:32,48:32}']=112*L(IDLE)
pads=Counter({r:n for r,n in u1['pad_hist'].items() if r})
comp['padding']=L(pads)
tot=sum(comp.values())
for k,v in comp.items(): print('  L %-62s %10.0f  (%.1f%%)'%(k,v,100*v/tot))
print('  total %.0f ; first-order a=2560/L=%.6e'%(tot,2560/tot))
