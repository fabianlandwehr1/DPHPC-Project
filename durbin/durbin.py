import numpy as np


def initialize(N, datatype=np.float32):
    r = np.fromfunction(lambda i: N + 1 - i, (N, ), dtype=datatype)
    return r


def kernel(r):

    y = np.empty_like(r)
    alpha = -r[0]
    beta = 1.0
    y[0] = -r[0]

    for k in range(1, r.shape[0]):
        beta *= 1.0 - alpha * alpha
        alpha = -(r[k] + np.dot(np.flip(r[:k]), y[:k])) / beta
        print(beta, alpha, np.dot(np.flip(r[:k]), y[:k]))
        y[:k] += alpha * np.flip(y[:k])
        y[k] = alpha

    return y


print(kernel(initialize(1024)))
