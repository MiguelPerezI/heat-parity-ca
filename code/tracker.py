"""Structure tracker for the 1-D rule u' = frac(mean3(u)+a), rounded. Finds active clusters per step,
links them in time by overlap, returns tracks with velocity and width."""
import numpy as np
DEC, TH = 9, 0.05
def step(u, a, bc):
    L, R = np.roll(u, 1), np.roll(u, -1)
    if bc == 'mirror': L[0] = u[1]; R[-1] = u[-2]
    return np.round(np.mod((L + u + R)/3.0 + a, 1.0), DEC)
def active(u, bc):
    d = np.mod(u - np.roll(u, -1), 1.0); d = np.minimum(d, 1-d)
    dl = np.roll(d, 1)
    if bc == 'mirror': d[-1] = 0; dl[0] = 0
    return np.maximum(d, dl) > TH
def clusters(A, gap=2):
    """list of (start, end) inclusive, merging gaps <= gap (no wrap handling: periodic runs are long enough)"""
    idx = np.flatnonzero(A)
    if idx.size == 0: return []
    out = []; s = p = idx[0]
    for i in idx[1:]:
        if i - p > gap + 1: out.append((s, p)); s = i
        p = i
    out.append((s, p)); return out
def run_tracks(u, a, bc, T, t0=0, maxw=40):
    """returns list of tracks: dict(t=[...], c=[...centre], w=[...width])"""
    tracks = []; live = []   # live: (track, last_centre)
    N = len(u)
    for t in range(T):
        if t >= t0:
            cl = [c for c in clusters(active(u, bc)) if c[1]-c[0]+1 <= maxw]
            cents = [0.5*(s+e) for s, e in cl]; widths = [e-s+1 for s, e in cl]
            newlive = []; used = set()
            for tr, lc in live:
                # nearest cluster within 3 cells (periodic distance for torus)
                best = None; bd = 3.1
                for j, c in enumerate(cents):
                    d = abs(c - lc); d = min(d, N - d) if bc == 'periodic' else d
                    if d < bd and j not in used: best, bd = j, d
                if best is not None:
                    used.add(best); tr['t'].append(t); tr['c'].append(cents[best]); tr['w'].append(widths[best]); newlive.append((tr, cents[best]))
            for j in range(len(cents)):
                if j not in used:
                    tr = {'t':[t], 'c':[cents[j]], 'w':[widths[j]]}; tracks.append(tr); newlive.append((tr, cents[j]))
            live = newlive
        u = step(u, a, bc)
    return tracks, u
def velocity(tr, N, bc):
    t = np.array(tr['t'], float); c = np.array(tr['c'], float)
    if bc == 'periodic': c = np.unwrap(c*2*np.pi/N)*N/(2*np.pi)   # unwrap periodic positions
    if len(t) < 2: return 0.0
    return np.polyfit(t, c, 1)[0]
