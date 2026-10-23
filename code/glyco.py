"""
Parity display of a Turing instability: the glycolysis reaction-diffusion model of
Cisneros-Ake, Gonzalez-Rodriguez, Gonzalez-Ramirez (Bull. Math. Biol. 84:100, 2022), Fig. 4:
  u_t = Du lap u - u + v(a+u^2),   v_t = Dv lap v + b - v(a+u^2),   Neumann on [0,1]^2,
  a=0.1, b=0.85, Du=0.0025, Dv=8 Du, dx=0.01 (101 nodes), dt=0.001, IC = equilibrium + 0.01*U[0,1].
Displayed as floor(q u) mod 2 with q=1e4. Tracks the quantities of Theorem 3.1:
q TV_h, E_theta P_h (Monte-Carlo and the exact formula), the aliasing deficit q TV_h - E_theta P_h,
the largest jump in quanta, and the 2x2 block entropy. Fits the growth rate of q TV_h in the linear
regime and the time t* at which the largest jump crosses one quantum.
usage: python3.12 glyco.py [T=3000]   (about 17 min for T=3000)
"""
import sys, os, numpy as np
from pngio import write_gray
here = os.path.dirname(os.path.abspath(__file__)); fig = here + "/../figures"

a, b, Du = 0.1, 0.85, 0.0025; Dv = 8*Du
N = 101; dx = 1.0/(N-1); dt = 1e-3; q = 1e4
T = float(sys.argv[1]) if len(sys.argv) > 1 else 3000
us, vs = b, b/(a+b*b)
rng = np.random.default_rng(1)
u = us + 0.01*rng.uniform(size=(N, N)); v = vs + 0.01*rng.uniform(size=(N, N))

def lap(w):
    P = np.pad(w, 1, mode='reflect')            # node-centred Neumann, ghost w_{-1} = w_1
    return (P[2:, 1:-1]+P[:-2, 1:-1]+P[1:-1, 2:]+P[1:-1, :-2]-4*w)/dx**2
def jumps(w): return np.concatenate([np.abs(np.diff(w, axis=0)).ravel(), np.abs(np.diff(w, axis=1)).ravel()])
def tri(x): y = np.mod(x, 2.0); return np.minimum(y, 2.0-y)
def parity(w, th=0.0): return np.floor(q*w+th).astype(np.int64) % 2
def perim(w, th):
    c = parity(w, th); return int(np.sum(c[1:] != c[:-1]) + np.sum(c[:, 1:] != c[:, :-1]))
def h2x2(w):
    c = parity(w); k = c[:-1, :-1]*8 + c[1:, :-1]*4 + c[:-1, 1:]*2 + c[1:, 1:]
    p = np.bincount(k.ravel(), minlength=16)/k.size; p = p[p > 0]
    return float(-(p*np.log2(p)).sum())
def mode(w):                                     # dominant Neumann cosine mode (DCT-I via even extension)
    E = np.concatenate([w, w[-2:0:-1]], 0); E = np.concatenate([E, E[:, -2:0:-1]], 1)
    F = np.abs(np.fft.fft2(E - E.mean()))[:N, :N]
    return tuple(int(i) for i in np.unravel_index(np.argmax(F), F.shape))

# linear theory on the grid: growth rate of cos(n pi x) cos(m pi y)
fu, fv = -1+2*b*b/(a+b*b), a+b*b; gu, gv = -2*b*b/(a+b*b), -(a+b*b)
def growth(n, m):
    k2 = (4/dx**2)*(np.sin(n*np.pi*dx/2)**2 + np.sin(m*np.pi*dx/2)**2)
    return np.linalg.eigvals(np.array([[fu-Du*k2, fv], [gu, gv-Dv*k2]])).real.max()
lam = {(n, m): growth(n, m) for n in range(10) for m in range(10)}
print(f"q={q:g}, |E|={2*N*(N-1)}; unstable modes (n,m): growth rate")
for k_, l_ in sorted(lam.items(), key=lambda kv: -kv[1]):
    if l_ > 0: print(f"  {k_}: {l_:.5f}")

def save_frames(w, t):
    s = 4                                        # 4x4 pixels per cell
    par = np.kron((255*(1-parity(w))).astype(np.uint8), np.ones((s, s), np.uint8))
    lo, hi = w.min(), w.max()                    # field rescaled to its own range, 32 grey levels
    g = np.floor(31.999*(w-lo)/(hi-lo)).astype(np.uint8)*8
    write_gray(f"{fig}/glyco_parity_t{t:04d}.png", par)
    write_gray(f"{fig}/glyco_field_t{t:04d}.png", np.kron(g, np.ones((s, s), np.uint8)))

table_t = [0, 10, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 800, 1000, 1500, 2000, 2500, 3000]
frame_t = {0, 200, 350, 1000}
rng2 = np.random.default_rng(0)
print(f"{'t':>6} {'range u':>9} {'bands':>6} {'mode':>7} {'qTV_h':>9} {'P(0)':>6} {'E_th P':>7} {'formula':>8} {'deficit':>8} {'max qJ':>7} {'H2x2':>5} {'flip MC':>8} {'flip fml':>8}")
dense = []                                       # (t, q TV_h, max q jump) every 5 time units
uprev = None                                     # state at the previous table time, for the flicker identity
steps = int(round(T/dt)); every = int(round(5/dt))
for k in range(steps+1):
    t = k*dt
    if k % every == 0:
        J = jumps(u); dense.append((t, q*J.sum(), q*J.max()))
        ti = int(round(t))
        if ti in table_t:
            f = tri(q*J).sum(); Pm = np.mean([perim(u, th) for th in rng2.uniform(0, 1, 64)])
            # flicker: fraction of cells whose colour changed since the previous table time (mean over phases),
            # against the time-edge version of the random-phase identity, mean_x phi(q |u(t,x) - u(s,x)|)
            if uprev is None: fl = flf = float('nan')
            else:
                fl = np.mean([(parity(u, th) != parity(uprev, th)).mean() for th in rng2.uniform(0, 1, 64)])
                flf = tri(q*np.abs(u-uprev)).mean()
            uprev = u.copy()
            print(f"{ti:6d} {u.max()-u.min():9.2e} {int(np.floor(q*u.max()))-int(np.floor(q*u.min())):6d} {str(mode(u)):>7} "
                  f"{q*J.sum():9.0f} {perim(u,0):6d} {Pm:7.0f} {f:7.0f} {q*J.sum()-f:8.0f} {q*J.max():7.2f} {h2x2(u):5.2f} {fl:8.4f} {flf:8.4f}", flush=True)
        if ti in frame_t: save_frames(u, ti)
    fuv = v*(a+u*u)
    u, v = u + dt*(Du*lap(u) - u + fuv), v + dt*(Dv*lap(v) + b - fuv)

D = np.array(dense)
w = (D[:, 0] >= 200) & (D[:, 0] <= 350)          # linear regime: stable modes gone, amplitude still small
slope = np.polyfit(D[w, 0], np.log(D[w, 1]), 1)[0]
i = np.argmax((D[:, 0] >= 100) & (D[:, 2] > 1))  # first crossing after the initial noise has diffused
t0, t1 = D[i-1, 0], D[i, 0]
tstar = t0 + (t1-t0)*np.log(1/D[i-1, 2])/np.log(D[i, 2]/D[i-1, 2])
print(f"growth rate of q TV_h on [200,350]: {slope:.5f}   (linear theory {lam[(3, 2)]:.5f})")
print(f"t* (max jump = one quantum): {tstar:.1f}")
print(f"final max jump {D[-1, 2]/q:.4e}: the display is aliased for every q > {q/D[-1, 2]:.1f}")
