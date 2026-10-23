"""Reflection principle: Neumann (mirror-ghost) box CA == periodic CA on the doubled domain
restricted to reflection-symmetric configurations. Checked exactly, 1-D and 2-D."""
import numpy as np
DEC = 9
def step1d(u, a, bc):
    L, R = np.roll(u, 1), np.roll(u, -1)
    if bc == 'neumann': L[0] = u[1]; R[-1] = u[-2]
    return np.round(np.mod((L + u + R)/3.0 + a, 1.0), DEC)
def double(u):  # even reflection about the wall nodes 0 and N-1 -> periodic of length 2(N-1)
    return np.concatenate([u, u[-2:0:-1]])
rng = np.random.default_rng(0); N = 400; a = 0.398
u = np.round(rng.random(N), DEC); v = double(u)
md = 0.0
for t in range(3000):
    u = step1d(u, a, 'neumann'); v = step1d(v, a, 'periodic')
    md = max(md, np.abs(v[:N] - u).max())
print(f"1-D: max |box - doubled periodic| over 3000 steps = {md:.3e}")

# 2-D, 5-point mean, mirror ghosts per axis vs doubled torus
def step2d(u, a, bc):
    def sh(ax, s):
        w = np.roll(u, s, ax)
        if bc == 'neumann':
            idx = [slice(None)]*2
            if s == 1:  idx[ax] = 0;  src = [slice(None)]*2; src[ax] = 1;  w[tuple(idx)] = u[tuple(src)]
            else:       idx[ax] = -1; src = [slice(None)]*2; src[ax] = -2; w[tuple(idx)] = u[tuple(src)]
        return w
    m = (u + sh(0,1) + sh(0,-1) + sh(1,1) + sh(1,-1))/5.0
    return np.round(np.mod(m + a, 1.0), DEC)
N = 64; u = np.round(rng.random((N,N)), DEC)
v = np.concatenate([u, u[-2:0:-1, :]], axis=0); v = np.concatenate([v, v[:, -2:0:-1]], axis=1)
md = 0.0
for t in range(1000):
    u = step2d(u, a, 'neumann'); v = step2d(v, a, 'periodic')
    md = max(md, np.abs(v[:N, :N] - u).max())
print(f"2-D: max |box - doubled torus| over 1000 steps = {md:.3e}   (doubled size {v.shape})")
