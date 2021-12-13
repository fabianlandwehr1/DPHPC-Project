import numpy as np

N = 1024
npt = 10


def initialize(N):
    data, radius = np.empty(N, dtype=np.float64), np.empty(N, dtype=np.float64)
    seed = 1
    for i in range(N):
        seed = (seed * 0x510f506f + 0xf59e9d39) & 0xFFFFFFFF
        data[i] = seed / 0xFFFFFFFF
        seed = (seed * 0x510f506f + 0xf59e9d39) & 0xFFFFFFFF
        radius[i] = seed / 0xFFFFFFFF
    return data, radius


def kernel(data, radius, npt):
    rmax = radius.max()
    res = np.zeros(npt, dtype=np.float64)
    for i in range(npt):
        r1 = rmax * i / npt
        r2 = rmax * (i + 1) / npt
        mask_r12 = np.logical_and((r1 <= radius), (radius < r2))
        values_r12 = data[mask_r12]
        res[i] = values_r12.mean()
    return res


d, r = initialize(N)
print(kernel(d, r, npt))
