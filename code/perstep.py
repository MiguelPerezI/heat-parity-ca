"""Per-step statistics in the bulk (cells 50..N-51), no tracking: structure density and wide-cluster density,
torus vs mirror box, 10 seeds x 10000 steps after 3000 transient."""
import numpy as np
from tracker import step, active, clusters
a, N, SEEDS, T, T0 = 0.398, 400, 10, 13000, 3000
for bc in ('periodic', 'mirror'):
    n_cl = n_wide = n_steps = 0; wall_wide = 0; wall_cl = 0
    for seed in range(SEEDS):
        rng = np.random.default_rng(100+seed); u = np.round(rng.random(N), 9)
        for t in range(T):
            if t >= T0:
                cl = clusters(active(u, bc)); n_steps += 1
                for s, e in cl:
                    c = 0.5*(s+e); w = e-s+1
                    if 50 <= c <= N-51:
                        n_cl += 1; n_wide += (w >= 6)
                    elif min(c, N-1-c) < 10:
                        wall_cl += 1; wall_wide += (w >= 6)
            u = step(u, a, bc)
    L = N-100
    print(f"{bc:8s}: bulk clusters/cell/step={n_cl/(n_steps*L):.4f}  wide(>=6)/cell/step={n_wide/(n_steps*L):.5f}  wide fraction={n_wide/n_cl:.4f}  "
          f"wide per cluster^2 (x1e3)={1e3*n_wide/(n_steps*L)/(n_cl/(n_steps*L))**2:.2f}  | within 10 cells of 0/N-1: clusters/cell/step={wall_cl/(n_steps*20):.4f} wide fraction={wall_wide/max(1,wall_cl):.4f}")
