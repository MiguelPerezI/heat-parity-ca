import numpy as np
from tracker import run_tracks, velocity
a, N = 0.398, 600
rng = np.random.default_rng(0); u = np.round(rng.random(N), 9)
tracks, _ = run_tracks(u, a, 'periodic', T=4000, t0=1000)
long = [tr for tr in tracks if len(tr['t']) >= 200]
print(f"tracks total {len(tracks)}, with lifetime >= 200 steps: {len(long)}")
vs = np.array([velocity(tr, N, 'periodic') for tr in long]); ws = np.array([np.mean(tr['w']) for tr in long]); ls = np.array([len(tr['t']) for tr in long])
hist, edges = np.histogram(vs, bins=np.arange(-1.05, 1.1, 0.1))
print("velocity histogram (cells/step):"); 
for h, e in zip(hist, edges[:-1]):
    if h: print(f"  [{e:+.2f},{e+0.1:+.2f}) : {h}")
print("velocity, mean width, lifetime of the 25 longest-lived:")
for i in np.argsort(-ls)[:25]: print(f"  v={vs[i]:+.3f}  w={ws[i]:5.1f}  life={ls[i]}")
