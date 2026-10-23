import sys; sys.path.insert(0,'.')
from analyze import run, out
import numpy as np
print("== late-time resolution test: t = 0.4, quantum 0.1 fixed, Jacobi Neumann ==")
for N in (64,128,256):
    dt = 2e-4*(64/N)**2; steps = int(round(0.4/dt))
    k,m,nb,e,H,c = run("jacobi", N=N, dt=dt, steps=steps, snaps=(steps,), tag="_late")[0]
    print(f"N={N:3d} steps={steps:6d} max={m:.3f} bands={nb} perimeter/N^2={e/N**2:.3f} H2x2={H:.3f} components={c}")
