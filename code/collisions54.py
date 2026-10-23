"""Rule 54 torus: outcomes of head-on w+ / w- collisions found by the tracker. A collision event: a w+ track and a w- track
end within 3 cells and 3 steps of each other. Outcome: what tracks start within 8 cells and <= 12 steps after:
0 new -> annihilation; only stationary -> gutter; two movers -> pass-through/exchange; otherwise mixed."""
import numpy as np
from eca import track, velocity, ether_words
ether_words(54, [1,0,0,0], 4)
N = 1400; out = {'annihilation':0, 'gutter':0, 'two movers':0, 'mixed':0}; n_ev = 0
def ic(rng):
    u = np.tile(np.array([1,0,0,0]), N//4 + 1)[:N]
    for _ in range(25):
        x = rng.integers(20, N-28); u[x:x+8] = rng.integers(0, 2, 8)
    return u
for seed in range(4):
    rng = np.random.default_rng(seed); u = ic(rng)
    tracks, _, _ = track(u, 54, 'periodic', 4, T=4000, t0=500, gap=2)
    T = [(tr, velocity(tr, N, 'periodic')) for tr in tracks if len(tr['t']) >= 20]
    plus = [(tr, v) for tr, v in T if v > 0.8]; minus = [(tr, v) for tr, v in T if v < -0.8]
    starts = [(tr['t'][0], tr['c'][0], ('mover' if abs(v) > 0.5 else 'stationary')) for tr, v in T]
    for trp, vp in plus:
        te, ce = trp['t'][-1], trp['c'][-1]
        for trm, vm in minus:
            if abs(trm['t'][-1] - te) <= 3 and abs(trm['c'][-1] - ce) <= 3:
                n_ev += 1
                new = [k for (t0, c0, k) in starts if te < t0 <= te + 12 and abs(c0 - ce) <= 8]
                if not new: out['annihilation'] += 1
                elif all(k == 'stationary' for k in new): out['gutter'] += 1
                elif new.count('mover') >= 2: out['two movers'] += 1
                else: out['mixed'] += 1
print(f"head-on w+/w- collisions found: {n_ev}; outcomes: {out}")
