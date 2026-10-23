import sys, numpy as np
from eca import track, velocity
rule, p, gap, N = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(0); u = rng.integers(0, 2, N)
tracks, u, dens = track(u, rule, 'periodic', p, T=4000, t0=1500, gap=gap)
sel = [tr for tr in tracks if len(tr['t']) >= 60]
vs = np.array([velocity(tr, N, 'periodic') for tr in sel]); ws = np.array([np.mean(tr['w']) for tr in sel]); ls = np.array([len(tr['t']) for tr in sel])
print(f"rule {rule}: defect density late = {dens[-500:].mean():.3f}; tracks >= 60 steps: {len(sel)}")
bins = np.arange(-1.0, 1.01, 0.05)
hist, edges = np.histogram(vs, bins=bins)
for h, e in zip(hist, edges[:-1]):
    if h: m = (vs>=e)&(vs<e+0.05); print(f"  v in [{e:+.2f},{e+0.05:+.2f}): {h:4d}  width {ws[m].mean():5.1f}  mean life {ls[m].mean():6.0f}  max life {ls[m].max()}")
