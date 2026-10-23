"""Rule 110: torus vs mirror box. Species composition (bulk-only tracks), per-step bulk defect/cluster density,
and defect density near the walls (is the wall a glider source?)."""
import sys, numpy as np
from eca import track, velocity, defects, clusters, step
rule, p, gap, N, T, T0, SEEDS = 110, 14, 3, 1400, 5000, 1500, 3
species = [(-0.5,'B'), (-0.25,'E'), (-0.111,'F'), (0.0,'C'), (0.2,'D'), (0.667,'A')]
def classify(v):
    k = int(np.argmin([abs(v - s) for s, _ in species])); return species[k][1] if abs(v - species[k][0]) < 0.05 else '?'
for bc in ('periodic', 'mirror'):
    comp = {}; nbulk = 0; dens_bulk = []; dens_wall = []; cl_bulk = []
    for seed in range(SEEDS):
        rng = np.random.default_rng(seed); u = rng.integers(0, 2, N)
        tracks, u, dens = track(u, rule, bc, p, T=T, t0=T0, gap=gap)
        for tr in tracks:
            if len(tr['t']) < 60: continue
            cm = np.mean(tr['c'])
            if not (150 <= cm <= N-150): continue
            s = classify(velocity(tr, N, bc)); comp[s] = comp.get(s, 0) + 1; nbulk += 1
        # per-step densities over the last 500 steps (continue the run)
        for t in range(500):
            D = defects(u, p, bc)
            dens_bulk.append(D[150:N-150].mean()); dens_wall.append(np.concatenate([D[:30], D[N-30:]]).mean())
            cl_bulk.append(len([c for c in clusters(D, gap) if 150 <= 0.5*(c[0]+c[1]) <= N-150]))
            u = step(u, rule, bc)
    tot = sum(comp.values())
    print(f"{bc:8s}: bulk tracks {tot}: " + ', '.join(f"{k}={v} ({v/tot:.2f})" for k, v in sorted(comp.items(), key=lambda kv: -kv[1])))
    print(f"          per-step: bulk defect density {np.mean(dens_bulk):.4f}, bulk clusters per 1100 cells {np.mean(cl_bulk):.1f}, defect density within 30 cells of the ends {np.mean(dens_wall):.4f}")
