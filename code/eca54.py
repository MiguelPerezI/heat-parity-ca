"""Rule 54 (left-right symmetric): sparse launch on its ether (tile 1000, temporal period 4), torus vs mirror box.
Box length 4k+1 so that both wall nodes sit on a '1' and the ether is mirror-compatible."""
import numpy as np
from eca import track, velocity, ether_words, defects, clusters, step
ether_words(54, [1,0,0,0], 4)
rule, p, gap, T, T0, SEEDS = 54, 4, 2, 5000, 1500, 3
def ic(N, rng, npatch=25):
    u = np.tile(np.array([1,0,0,0]), N//4 + 1)[:N]
    for _ in range(npatch):
        x = rng.integers(20, N-28); u[x:x+8] = rng.integers(0, 2, 8)
    return u
for bc, N in (('periodic', 1400), ('mirror', 1401)):
    vs_all = []; ws_all = []; ls_all = []; dens_bulk = []; dens_wall = []; cl_bulk = []
    for seed in range(SEEDS):
        rng = np.random.default_rng(seed); u = ic(N, rng)
        tracks, u, dens = track(u, rule, bc, p, T=T, t0=T0, gap=gap)
        for tr in tracks:
            if len(tr['t']) < 60: continue
            cm = np.mean(tr['c'])
            if not (150 <= cm <= N-150): continue
            vs_all.append(velocity(tr, N, bc)); ws_all.append(np.mean(tr['w'])); ls_all.append(len(tr['t']))
        for t in range(500):
            D = defects(u, p, bc, rule)
            dens_bulk.append(D[150:N-150].mean()); dens_wall.append(np.concatenate([D[:30], D[N-30:]]).mean())
            cl_bulk.append(len([c for c in clusters(D, gap) if 150 <= 0.5*(c[0]+c[1]) <= N-150]))
            u = step(u, rule, bc)
    vs = np.array(vs_all); ws = np.array(ws_all); ls = np.array(ls_all)
    print(f"{bc:8s}: bulk tracks {len(vs)}; per-step bulk defect density {np.mean(dens_bulk):.4f}, bulk clusters {np.mean(cl_bulk):.1f}, defect density within 30 cells of ends {np.mean(dens_wall):.4f}")
    hist, edges = np.histogram(vs, bins=np.arange(-1.05, 1.1, 0.1))
    for h, e in zip(hist, edges[:-1]):
        if h: m = (vs>=e)&(vs<e+0.1); print(f"    v in [{e:+.2f},{e+0.1:+.2f}): {h:4d} ({h/len(vs):.2f})  width {ws[m].mean():5.1f}  mean life {ls[m].mean():5.0f}")
