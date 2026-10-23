"""Elementary cellular automata on torus / mirror box, periodic-background (ether) detection, glider tracking."""
import numpy as np
def step(u, rule, bc):
    L, R = np.roll(u, 1), np.roll(u, -1)
    if bc == 'mirror': L[0] = u[1]; R[-1] = u[-2]
    return (rule >> (4*L + 2*u + R)) & 1
def words(u, p):
    """integer code of the p-word starting at each cell (periodic wrap for the code only)"""
    w = np.zeros(len(u), np.int64)
    for k in range(p): w = 2*w + np.roll(u, -k)
    return w
def ether_set(u, p):
    """most frequent p-word and its p cyclic rotations (the background at this time)"""
    w = words(u, p); vals, cnt = np.unique(w, return_counts=True); best = vals[np.argmax(cnt)]
    bits = [(best >> (p-1-k)) & 1 for k in range(p)]
    rots = set()
    for r in range(p):
        b = bits[r:] + bits[:r]; v = 0
        for x in b: v = 2*v + x
        rots.add(v)
    return rots
FIXED = {}   # rule -> set of all ether words (all temporal phases); filled by ether_words()
def ether_words(rule, tile, period):
    """all p-words of the ether over one temporal period, from a pure tiling of `tile`"""
    p = len(tile); u = np.tile(np.array(tile), 20); S = set()
    for t in range(period):
        S |= set(np.unique(words(u, p)).tolist()); u = step(u, rule, 'periodic')
    assert np.array_equal(u, np.tile(np.array(tile), 20)), "tile is not periodic with that period"
    FIXED[rule] = S; return S
def defects(u, p, bc, rule=None):
    """cells covered by a non-ether p-word, flagged symmetrically (word centred on the cell).
    Uses FIXED[rule] if available, else the rotations of the most frequent word at this time."""
    S = FIXED.get(rule) if rule in FIXED else ether_set(u, p); w = words(u, p)
    D0 = np.array([x not in S for x in w])          # word starting at x is non-ether
    N = len(u); D = np.roll(D0, p//2)                 # centre the flag
    if bc == 'mirror':
        D[:p//2] = D0[0]                              # ends: no wrapped information; copy the nearest valid flag
        D[N-p+1+p//2:] = D0[N-p]
    return D
def clusters(D, gap):
    idx = np.flatnonzero(D)
    if idx.size == 0: return []
    out = []; s = p = idx[0]
    for i in idx[1:]:
        if i - p > gap + 1: out.append((s, p)); s = i
        p = i
    out.append((s, p)); return out
def track(u, rule, bc, p, T, t0, gap, maxw=60, match=4):
    rule_ = rule
    N = len(u); tracks = []; live = []; dens = []
    for t in range(T):
        if t >= t0:
            D = defects(u, p, bc, rule_); dens.append(D.mean())
            cl = [c for c in clusters(D, gap) if c[1]-c[0]+1 <= maxw]
            cents = [0.5*(s+e) for s, e in cl]; wid = [e-s+1 for s, e in cl]
            newlive = []; used = set()
            for tr, lc in live:
                best = None; bd = match + 0.1
                for j, c in enumerate(cents):
                    d = abs(c-lc); d = min(d, N-d) if bc == 'periodic' else d
                    if d < bd and j not in used: best, bd = j, d
                if best is not None:
                    used.add(best); tr['t'].append(t); tr['c'].append(cents[best]); tr['w'].append(wid[best]); newlive.append((tr, cents[best]))
            for j in range(len(cents)):
                if j not in used:
                    tr = {'t':[t], 'c':[cents[j]], 'w':[wid[j]]}; tracks.append(tr); newlive.append((tr, cents[j]))
            live = newlive
        u = step(u, rule, bc)
    return tracks, u, np.array(dens)
def velocity(tr, N, bc):
    t = np.array(tr['t'], float); c = np.array(tr['c'], float)
    if bc == 'periodic': c = np.unwrap(c*2*np.pi/N)*N/(2*np.pi)
    return np.polyfit(t, c, 1)[0] if len(t) > 1 else 0.0
