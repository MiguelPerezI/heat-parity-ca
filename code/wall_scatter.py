"""Rule 110 mirror box: what goes into each wall and what comes out. Tracks (>=30 steps) that end within 40 cells of a wall
moving toward it = incoming; tracks that start within 40 cells moving away = outgoing. Also total in/out per wall."""
import sys, numpy as np
from eca import track, velocity, ether_words
ether_words(54, [1,0,0,0], 4); ether_words(110, [1,1,0,0,0,1,0,0,1,1,0,1,1,1], 7)
rule, p, gap = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
N, T, T0, SEEDS, W = 1400 + (1 if rule == 54 else 0), 6000, 1500, 3, 40
species = {110: [(-0.5,'B'), (-0.25,'E'), (-0.111,'F'), (0.0,'C'), (0.2,'D'), (0.667,'A')],
            54: [(-1.0,'w-'), (-0.5,'?-'), (0.0,'g'), (0.5,'?+'), (1.0,'w+')]}[rule]
def classify(v):
    k = int(np.argmin([abs(v - s) for s, _ in species])); return species[k][1] if abs(v - species[k][0]) < 0.07 else '?'
inc = {'L':{}, 'R':{}}; out = {'L':{}, 'R':{}}
def ic(rng):
    if rule == 110: return rng.integers(0, 2, N)
    u = np.tile(np.array([1,0,0,0]), N//4 + 1)[:N]
    for _ in range(25):
        x = rng.integers(20, N-28); u[x:x+8] = rng.integers(0, 2, 8)
    return u
for seed in range(SEEDS):
    rng = np.random.default_rng(seed); u = ic(rng)
    tracks, u, dens = track(u, rule, 'mirror', p, T=T, t0=T0, gap=gap)
    for tr in tracks:
        if len(tr['t']) < 30: continue
        v = velocity(tr, N, 'mirror'); s = classify(v); c0, c1 = tr['c'][0], tr['c'][-1]
        if abs(v) < 0.03: continue
        if c1 < W and v < 0: inc['L'][s] = inc['L'].get(s, 0) + 1
        if c1 > N-1-W and v > 0: inc['R'][s] = inc['R'].get(s, 0) + 1
        if c0 < W and v > 0: out['L'][s] = out['L'].get(s, 0) + 1
        if c0 > N-1-W and v < 0: out['R'][s] = out['R'].get(s, 0) + 1
for side in ('L', 'R'):
    print(f"rule {rule} wall {side}: incoming {sum(inc[side].values()):4d} {dict(sorted(inc[side].items()))}")
    print(f"rule {rule} wall {side}: outgoing {sum(out[side].values()):4d} {dict(sorted(out[side].items()))}")
