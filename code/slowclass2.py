import numpy as np
from tracker import run_tracks, velocity
a, N, SEEDS = 0.398, 400, 5
for bc in ('periodic', 'mirror'):
    bulk = {'fast':0, 'slow':0}; wallhist = np.zeros(4); allslow = 0; slow_life = []; fast_life = []
    for seed in range(SEEDS):
        rng = np.random.default_rng(seed); u = np.round(rng.random(N), 9)
        tracks, _ = run_tracks(u, a, bc, T=10000, t0=3000)
        for tr in tracks:
            if len(tr['t']) < 50: continue
            v = abs(velocity(tr, N, bc)); cmean = np.mean(tr['c'])
            d = min(cmean, N-1-cmean)
            cls = 'slow' if v < 0.05 else 'fast'
            (slow_life if cls=='slow' else fast_life).append(len(tr['t']))
            if 50 <= cmean <= N-50: bulk[cls] += 1           # same bulk criterion for both geometries
            if cls == 'slow':
                allslow += 1; wallhist[0 if d < 5 else 1 if d < 20 else 2 if d < 50 else 3] += 1
    n = bulk['fast'] + bulk['slow']
    print(f"{bc:8s}: bulk tracks {n}, slow fraction in bulk = {bulk['slow']/n:.3f} ({bulk['slow']} slow, {bulk['fast']} fast); "
          f"mean life slow {np.mean(slow_life):.0f}, fast {np.mean(fast_life):.0f}")
    cells = np.array([5, 15, 30, N/2-50])*2   # cells per distance bin (both sides)
    print(f"          slow-track centres by distance to nearest wall [0-5,5-20,20-50,>50]: counts {wallhist.astype(int)}, per cell {np.round(wallhist/cells,3)}")
