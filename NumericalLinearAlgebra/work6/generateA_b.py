import numpy as np

rng = np.random.default_rng(seed=0)
N=[10, 20, 40, 80]
A = [[] for _ in range(len(N))]
for n in N:
    for i in range(100):
        _A = rng.random((n, n))
        A[N.index(n)].append(_A)

# A= A[N.index(n)][i]となる