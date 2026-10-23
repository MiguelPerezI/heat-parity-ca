"""
Continuous cellular automaton 'diffusion mod 1' (Wolfram, NKS ch. 4):
    u' = frac( u + sum_m r_m * d2_m u + a )
on a box in any dimension, with boundary handled by the ternary-cube rule:
  periodic   : wrap
  neumann    : reflected ghost (homogeneous Neumann), per axis, per trit
With r_m = 1/(2n+1) this is 'mean of the (2n+1)-point neighbourhood plus a, mod 1'.
"""
import numpy as np

def d2_axis(u, m, bc):
    if bc == 'periodic':
        return np.roll(u, -1, axis=m) - 2*u + np.roll(u, 1, axis=m)
    # neumann: ghost = inner neighbour (trit 0: u_{-1}=u_{1}; trit 2: u_{N}=u_{N-2})
    up = np.roll(u, -1, axis=m); dn = np.roll(u, 1, axis=m)
    sl_first = [slice(None)]*u.ndim; sl_first[m] = 0
    sl_last  = [slice(None)]*u.ndim; sl_last[m]  = -1
    sl_1 = [slice(None)]*u.ndim; sl_1[m] = 1
    sl_m2 = [slice(None)]*u.ndim; sl_m2[m] = -2
    dn[tuple(sl_first)] = u[tuple(sl_1)]
    up[tuple(sl_last)]  = u[tuple(sl_m2)]
    return up - 2*u + dn

def step(u, a, r, bc):
    lap = sum(r*d2_axis(u, m, bc) for m in range(u.ndim))
    return np.mod(u + lap + a, 1.0)

def active(u, bc, eps=1e-6):
    """cells whose neighbourhood is not flat (Laplacian != 0)"""
    lap = sum(d2_axis(u, m, bc) for m in range(u.ndim))
    return np.abs(lap) > eps

def trits(shape):
    """ternary word array: per axis 0 (lower wall), 1 (free), 2 (upper wall)"""
    grids = np.meshgrid(*[np.arange(N) for N in shape], indexing='ij')
    return [np.where(g == 0, 0, np.where(g == N-1, 2, 1)) for g, N in zip(grids, shape)]

def codim(shape):
    return sum((t != 1).astype(int) for t in trits(shape))
