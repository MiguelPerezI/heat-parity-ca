"""Rule 54, sparse launches with the exact mover counter: one random 8-cell patch in the centre of the ether,
torus of 400 cells vs mirror box of 401 cells, same patches. Movers N+ + N- and off-ether cells over time.
Patches that leave a stationary gun (off-ether cells growing) are reported separately."""
import numpy as np
from eca import step
from rule54_charge import counts, ether
times = (100, 300, 600, 1000, 1500, 2000)
res = {'periodic': [], 'mirror': []}; guns = {'periodic': 0, 'mirror': 0}; S = 200
for seed in range(S):
    patch = np.random.default_rng(seed).integers(0, 2, 8)
    for bc, N in (('periodic', 400), ('mirror', 401)):
        u = ether(N); u[196:204] = patch; E = ether(N); row = []
        for t in range(max(times) + 1):
            if t in times: row.append(counts(u, E, bc))
            u = step(u, 54, bc); E = step(E, 54, bc)
        res[bc].append(row)
R = {bc: np.array(v) for bc, v in res.items()}
np.savez("rule54_sparse.npz", **R)
for bc in R:
    mv = R[bc][:, :, 0] + R[bc][:, :, 1]; off = R[bc][:, :, 2]
    quiet = off[:, -1] == 0
    print(f"{bc:8s}: mean movers at t={times}: {np.round(mv.mean(axis=0), 3).tolist()}; "
          f"seeds with movers left at t=2000: {(mv[:, -1] > 0).sum()}/{S}; back to pure ether: {quiet.sum()}/{S}; "
          f"stationary debris only: {((mv[:, -1] == 0) & ~quiet).sum()}/{S}")
mvp = R['periodic'][:, -1, 0] + R['periodic'][:, -1, 1]
print("torus survivors at t=2000, histogram:", np.bincount(mvp).tolist())
imb = np.abs(R['periodic'][:, 0, 0] - R['periodic'][:, 0, 1])
print("torus: |N+ - N-| at t=100 vs movers at t=2000, agreement in", int((imb == mvp).sum()), "of", S, "seeds")

gun = (R['periodic'][:, -1, 0] + R['periodic'][:, -1, 1] > 10) | (R['mirror'][:, -1, 0] + R['mirror'][:, -1, 1] > 10)
ng = ~gun; print(f"\nwithout guns ({ng.sum()} seeds):")
for bc in R:
    mv = R[bc][ng, :, 0] + R[bc][ng, :, 1]
    print(f"  {bc:8s}: mean movers {np.round(mv.mean(axis=0), 3).tolist()}, seeds with movers at t=2000: {(mv[:, -1] > 0).sum()}")
P = R['periodic'][ng]
fin = P[:, -1, 0] + P[:, -1, 1]; d0 = P[:, 1, 0] - P[:, 1, 1]; dT = P[:, -1, 0] - P[:, -1, 1]
print("  torus survivors: all same direction in", int(((P[:, -1, 0] == 0) | (P[:, -1, 1] == 0)).sum()), "of", ng.sum(),
      "; N+ - N- at t=300 equals N+ - N- at t=2000 in", int((d0 == dT).sum()), "seeds")
