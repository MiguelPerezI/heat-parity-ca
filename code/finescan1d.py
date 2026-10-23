"""Fine scan of u' = frac(k*mean3(u) + a) in 1-D, periodic, for class-4 candidates."""
import numpy as np
from contca import active
N, T = 256, 400
def step_k(u, a, k):
    m = (np.roll(u,1)+u+np.roll(u,-1))/3.0
    return np.mod(k*m + a, 1.0)
rng = np.random.default_rng(2)
print(f"{'k':>4} {'a':>6} {'w200':>5} {'w400':>5} {'slope':>6} {'inner':>6} {'actR':>6} {'ncompR':>6} {'tag':>8}")
for k in (1.0, 1.5, 2.0):
    for a in np.round(np.arange(0.0, 1.0, 0.005), 3):
        u = np.zeros(N); u[N//2] = 0.5
        w = {}
        for t in range(T+1):
            if t in (200, 400):
                A = active(u, 'periodic'); w[t] = int(A.sum())
                if t == 400:
                    idx = np.flatnonzero(A); span = (idx.max()-idx.min()+1) if idx.size else 0
                    inner = A.sum()/span if span else 0.0
            u = step_k(u, a, k)
        slope = (w[400]-w[200])/200.0
        # random IC
        u = rng.random(N); acts=[]; nc=[]
        for t in range(T+1):
            if t >= 200:
                A = active(u,'periodic'); acts.append(A.mean())
                d = np.diff(np.concatenate([[0],A.astype(int),[0]])); nc.append(int((d==1).sum()))
            u = step_k(u, a, k)
        actR = float(np.mean(acts)); ncR = float(np.mean(nc))
        if w[400] == 0: tag = 'dies'
        elif slope < 0.05 and w[400] < 40: tag = 'LOCAL'
        elif slope < 0.2: tag = 'slow'
        else: tag = 'spread'
        if tag in ('LOCAL',) or (0.02 < actR < 0.6):
            print(f"{k:4.1f} {a:6.3f} {w[200]:5d} {w[400]:5d} {slope:6.2f} {inner:6.2f} {actR:6.3f} {ncR:6.1f} {tag:>8}")
