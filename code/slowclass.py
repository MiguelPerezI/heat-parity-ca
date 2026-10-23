"""Does the box change the composition of the structure population? Fraction of slow/wide tracks (|v|<0.05)
and velocity-class counts, torus vs mirror box, 3 seeds, late window."""
import numpy as np
from tracker import run_tracks, velocity
a, N = 0.398, 400
for bc in ('periodic', 'mirror'):
    tot = {'fast':0, 'slow':0}; lifes = []
    for seed in range(3):
        rng = np.random.default_rng(seed); u = np.round(rng.random(N), 9)
        tracks, _ = run_tracks(u, a, bc, T=8000, t0=3000)
        sel = [tr for tr in tracks if len(tr['t']) >= 50]
        for tr in sel:
            v = abs(velocity(tr, N, bc)); lifes.append(len(tr['t']))
            if bc == 'mirror' and (min(tr['c']) < 5 or max(tr['c']) > N-6): continue   # tracks touching the wall get bent; skip for class
            tot['fast' if v >= 0.05 else 'slow'] += 1
    n = tot['fast'] + tot['slow']
    print(f"{bc:8s}: tracks>=50 steps: {n}, slow/wide fraction = {tot['slow']/n:.3f} ({tot['slow']} slow, {tot['fast']} fast), mean lifetime {np.mean(lifes):.1f}, max {max(lifes)}")
