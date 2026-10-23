"""Rule 54: movers launched from one central patch; number alive vs time, torus (400) vs mirror box (401), 12 seeds.
Seeds whose patch nucleates a spreading phase (>20 'movers') are reported separately."""
import numpy as np
from eca import step, ether_words, track, velocity
ether_words(54, [1,0,0,0], 4)
def count_movers(box, bc, T=160):
    M = len(box); tr, _, _ = track(box.copy(), 54, bc, 4, T=T, t0=0, gap=2)
    return sum(1 for x in tr if len(x['t']) >= 80 and abs(velocity(x, M, bc)) > 0.5)
checks = (60, 300, 600, 900, 1500)
for bc, M in (('periodic', 400), ('mirror', 401)):
    normal = []; invaded = 0
    for seed in range(12):
        rng = np.random.default_rng(seed)
        box = np.tile(np.array([1,0,0,0]), M//4 + 1)[:M]; box[196:204] = rng.integers(0, 2, 8)
        c = {}
        for k in range(1501):
            if k in checks: c[k] = count_movers(box, bc)
            box = step(box, 54, bc)
        if c[60] > 20: invaded += 1; continue
        normal.append([c[k] for k in checks])
    A = np.array(normal)
    print(f"{bc:8s}: {len(normal)} normal seeds ({invaded} phase-invaded excluded); movers alive at t={checks}: mean {np.round(A.mean(axis=0),2)}; per seed {A.tolist()}")
