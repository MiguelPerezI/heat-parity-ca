"""
Analysis of the 'parity of floor(10 u)' pattern of main21D.py.
(a) reproduces main21D.py dynamics (in-place sweep, 3 frozen walls, buggy left wall),
(b) clean Jacobi scheme with homogeneous Neumann on all walls,
and tracks complexity measures of the binary pattern b = floor(q*u) mod 2.
"""
import numpy as np, zlib, struct, sys, os
out = os.path.dirname(os.path.abspath(__file__)) + "/../out"; os.makedirs(out, exist_ok=True)

def png(path, b, scale=4):
    a = (255*(1-b)).astype(np.uint8)              # white = parity 1 as in the script
    a = np.kron(a, np.ones((scale,scale),np.uint8))
    h, w = a.shape
    raw = b''.join(b'\x00'+a[y].tobytes() for y in range(h))
    def chunk(t,d): return struct.pack('>I',len(d))+t+d+struct.pack('>I',zlib.crc32(t+d)&0xffffffff)
    open(path,'wb').write(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',w,h,8,0,0,0,0))
                          +chunk(b'IDAT',zlib.compress(raw,9))+chunk(b'IEND',b''))

def measures(b):
    edges = np.sum(b[1:,:]!=b[:-1,:]) + np.sum(b[:,1:]!=b[:,:-1])   # perimeter (boundary edges)
    # 2x2 block entropy
    blocks = b[:-1,:-1]*8+b[:-1,1:]*4+b[1:,:-1]*2+b[1:,1:]
    p = np.bincount(blocks.ravel(), minlength=16)/blocks.size
    H = -np.sum(p[p>0]*np.log2(p[p>0]))
    # connected components (4-neighbour), union-find
    n = b.shape[0]*b.shape[1]; parent = np.arange(n)
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]]; x=parent[x]
        return x
    H_,W_ = b.shape
    for y in range(H_):
        for x in range(W_):
            i=y*W_+x
            if x+1<W_ and b[y,x]==b[y,x+1]:
                a,c=find(i),find(i+1)
                if a!=c: parent[a]=c
            if y+1<H_ and b[y,x]==b[y+1,x]:
                a,c=find(i),find(i+W_)
                if a!=c: parent[a]=c
    comps = len(set(find(i) for i in range(n)))
    return edges, H, comps

def run(mode, N=64, q=10.0, kd=0.2, dt=2e-4, L=1.0, steps=3000, snaps=(0,1,5,20,60,150,400,1000,3000), tag=""):
    dx = L/N; udx2 = 1/dx**2; r = kd*dt*udx2
    u = np.zeros((N,N)); s = slice(int(20*N/64), int(44*N/64)); u[s,s] = 100.0
    rows=[]
    def lap_jacobi(u):
        # homogeneous Neumann via ghost reflection on all walls (ternary-cube rule)
        up = np.empty_like(u)
        up[1:-1,1:-1] = u[2:,1:-1]+u[:-2,1:-1]+u[1:-1,2:]+u[1:-1,:-2]-4*u[1:-1,1:-1]
        # recompute with padded reflection for simplicity
        P = np.pad(u,1,mode='reflect')
        return P[2:,1:-1]+P[:-2,1:-1]+P[1:-1,2:]+P[1:-1,:-2]-4*u
    for k in range(steps+1):
        if k in snaps:
            b = (np.floor(q*u).astype(np.int64) % 2)
            e,H,c = measures(b); nb = int(np.floor(q*u.max()))-int(np.floor(q*u.min()))
            rows.append((k, u.max(), nb, e, H, c))
            png(f"{out}/{mode}{tag}_N{N}_k{k:05d}.png", b, scale=max(1,256//N))
        if mode=="script":
            # in-place sweep exactly as main21D.py after the first press (u is v):
            # interior updated row by row reading already-updated neighbours; walls frozen except left wall bug
            v = u  # alias
            for jj in range(1,N-1):
                for ii in range(1,N-1):
                    lap = (v[jj,ii-1]-2*v[jj,ii]+v[jj,ii+1])*udx2 + (v[jj-1,ii]-2*v[jj,ii]+v[jj+1,ii])*udx2
                    v[jj,ii] = v[jj,ii] + dt*kd*lap
            # left wall: w = u[43] + dt*k*lap  (global i == 43 bug) ; u[43] is cell (0,43)
            for ii in range(0,N-1):
                lap = (v[ii+1,0]-2*v[ii,0]+v[ii-1,0])*udx2 + (2*v[ii,1]-2*v[ii,0])*udx2
                v[ii,0] = v[0,43] + dt*kd*lap
        else:
            u = u + r*lap_jacobi(u)/1.0 * 1.0
            u = u  # Jacobi, Neumann all walls
    return rows

def table(rows, title):
    print(f"\n== {title} ==")
    print(f"{'step':>6} {'max u':>9} {'bands':>6} {'perimeter':>10} {'H2x2':>6} {'components':>11}")
    for k,m,nb,e,H,c in rows: print(f"{k:6d} {m:9.3f} {nb:6d} {e:10d} {H:6.3f} {c:11d}")

if __name__=="__main__":
    snaps=(0,1,5,20,60,150,400,1000,2000,3000)
    table(run("script", steps=3000, snaps=snaps), "main21D.py dynamics (in-place sweep, frozen walls), q=10, N=64")
    table(run("jacobi", steps=3000, snaps=snaps), "clean Jacobi, Neumann all walls, q=10, N=64")
    # resolution test at fixed physical time t = 150*2e-4 = 0.03: N=64,128,256 (dt scaled by (64/N)^2)
    print("\n== resolution test, clean Jacobi, same physical time t=0.03, same quantum 0.1 ==")
    for N in (64,128,256):
        dt = 2e-4*(64/N)**2; steps = int(round(0.03/dt))
        rows = run("jacobi", N=N, dt=dt, steps=steps, snaps=(steps,), tag="_res")
        k,m,nb,e,H,c = rows[0]
        print(f"N={N:3d}  perimeter/N^2={e/N**2:.3f}  H2x2={H:.3f}  components/N^2={c/N**2:.4f}  bands={nb}")
    print("\n== resolution test with quantum scaled as 0.1*(64/N) (aliasing ratio kept) ==")
    for N in (64,128,256):
        dt = 2e-4*(64/N)**2; steps = int(round(0.03/dt))
        rows = run("jacobi", N=N, q=10.0*N/64, dt=dt, steps=steps, snaps=(steps,), tag="_qscaled")
        k,m,nb,e,H,c = rows[0]
        print(f"N={N:3d}  perimeter/N^2={e/N**2:.3f}  H2x2={H:.3f}  components/N^2={c/N**2:.4f}  bands={nb}")
