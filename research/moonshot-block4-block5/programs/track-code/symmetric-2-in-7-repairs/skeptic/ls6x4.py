import sys; sys.path.insert(0,'')
import numpy as np
def res(a):
    A=np.array([[1,a,0,0],[a,1,0,0],[1,0,0,a],[a,0,0,1],[0,1,a,0],[0,a,1,0]],dtype=complex)
    t=np.zeros(6);t[0]=1
    b=np.linalg.lstsq(A,t,rcond=None)[0]; r=A@b-t
    return np.sum(abs(r)**2), np.round(abs(r)**2,4)
for a in [1,-1,0.5,0.9,1.1,0.0,2,1j, 0.7071]:
    print(a,res(a))
xs=np.linspace(-3,3,6001); v=[res(x)[0] for x in xs]; i=int(np.argmin(v)); print('min real a',xs[i],v[i])
