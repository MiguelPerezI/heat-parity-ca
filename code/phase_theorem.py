"""
Test of the random-phase identity for the parity perimeter:
  E_theta P_h = sum_edges dist(q*|du_e|, 2Z)      (triangle wave of the edge jump)
and of discrete TVD:  sum_edges |du_e| nonincreasing.
Jacobi scheme, Neumann (reflected ghosts), 64x64, same setup as main21D.py.
"""
import numpy as np
N=64; L=1.0; kd=0.2; dt=2e-4; q=10.0
dx=L/N; r=kd*dt/dx**2
u=np.zeros((N,N)); u[20:44,20:44]=100.0
rng=np.random.default_rng(0)
def lap(u):
    P=np.pad(u,1,mode='reflect'); return P[2:,1:-1]+P[:-2,1:-1]+P[1:-1,2:]+P[1:-1,:-2]-4*u
def jumps(u): return np.concatenate([np.abs(u[1:,:]-u[:-1,:]).ravel(), np.abs(u[:,1:]-u[:,:-1]).ravel()])
def perim(u,theta):
    b=np.floor(q*u+theta).astype(np.int64)%2
    return np.sum(b[1:,:]!=b[:-1,:])+np.sum(b[:,1:]!=b[:,:-1])
def tri(x):  # distance to nearest even integer
    y=np.mod(x,2.0); return np.minimum(y,2.0-y)
print(f"{'step':>6} {'q*TV_h':>10} {'P(theta=0)':>11} {'mean_theta P':>13} {'formula':>9} {'max q*jump':>11}")
checks=(0,1,5,20,60,150,400,1000,1500,2000,2500,3000,4000)
tv_prev=None
for k in range(4001):
    if k in checks:
        J=jumps(u); tv=J.sum()
        assert tv_prev is None or tv<=tv_prev+1e-9, "TVD violated"
        tv_prev=tv
        Pmean=np.mean([perim(u,th) for th in rng.uniform(0,1,200)])
        print(f"{k:6d} {q*tv:10.1f} {perim(u,0.0):11d} {Pmean:13.1f} {tri(q*J).sum():9.1f} {q*J.max():11.2f}")
    u=u+r*lap(u)
