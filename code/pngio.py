import numpy as np, zlib, struct
def write_gray(path, a):
    a = np.asarray(a, dtype=np.uint8); h, w = a.shape
    raw = b''.join(b'\x00' + a[y].tobytes() for y in range(h))
    def chunk(t, d): return struct.pack('>I', len(d)) + t + d + struct.pack('>I', zlib.crc32(t+d) & 0xffffffff)
    open(path, 'wb').write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 0, 0, 0, 0))
                           + chunk(b'IDAT', zlib.compress(raw, 9)) + chunk(b'IEND', b''))
def montage(tiles, cols=8, gap=4, bg=128):
    h, w = tiles[0].shape; rows = (len(tiles)+cols-1)//cols
    M = np.full((rows*(h+gap), cols*(w+gap)), bg, np.uint8)
    for i, t in enumerate(tiles):
        r, c = divmod(i, cols); M[r*(h+gap):r*(h+gap)+h, c*(w+gap):c*(w+gap)+w] = t
    return M
def upscale(a, s): return np.kron(a, np.ones((s, s), np.uint8))
