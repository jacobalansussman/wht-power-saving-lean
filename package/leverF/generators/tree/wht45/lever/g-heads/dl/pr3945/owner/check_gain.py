"""Exact arithmetic certificate for OWNER-LIFT-PROOF.txt's uniform gain."""
from fractions import Fraction as Q
import json
def need(ok,why):
    if not ok:raise ValueError(why)
def log_bounds(x,m):
    lower=2*sum((x**(2*j+1)/Q(2*j+1) for j in range(m)),Q(0))
    upper=lower+2*x**(2*m+1)/((2*m+1)*(1-x*x))
    return lower,upper
a,b=log_bounds(Q(1,3),3)
need(Q(693,1000)<a and b<Q(694,1000),'log 2 interval')
c,d=log_bounds(Q(1,2),6)
need(Q(1098,1000)<c and d<Q(1099,1000),'log 3 interval')
derivative_lower=7*Q(693,500)*Q(498612,125000)-26*Q(347,500)-18*Q(1099,1000)
need(derivative_lower==Q(13660453,15625000) and derivative_lower>0,'derivative lower bound')
need(16-13*2-6*3+7*4==0,'D(1)')
need(144-117*2-54*3+63*4==0,'rank mass')
need(9*derivative_lower==Q(122944077,15625000),'strict moment margin')
print(json.dumps({'status':'PASS_EXACT_UNIFORM_MOMENT_GAIN','kappa_domain':'0 < kappa <= 1/500','gain_per_replica_at_least':'(122944077/15625000)*kappa','derivative_lower':str(derivative_lower),'scope':'Rational log-series and derivative bounds; literal profile and physical admission are separately source checked.'}),flush=True)
