"""Structure count vs time: torus vs mirror box, same seeds. If walls remove structure types, the box count decays faster."""
import numpy as np
from tracker import step, active, clusters
a, N, T = 0.398, 400, 20000
checks = [0, 500, 1000, 2000, 4000, 8000, 12000, 16000, 19999]
for bc in ('periodic', 'mirror'):
    rows = []
    for seed in range(3):
        rng = np.random.default_rng(seed); u = np.round(rng.random(N), 9); counts = {}
        for t in range(T):
            if t in checks: counts[t] = len(clusters(active(u, bc)))
            u = step(u, a, bc)
        rows.append([counts[t] for t in checks])
    m = np.mean(rows, axis=0)
    print(f"{bc:8s} clusters at t={checks}:"); print("         ", ' '.join(f"{x:6.1f}" for x in m))
