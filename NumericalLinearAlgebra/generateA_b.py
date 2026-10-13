import numpy as np

rng = np.random.default_rng(seed=0)
N=[100, 200, 400, 800]
A_b = [[] for _ in range(len(N))]
for n in N:
    for i in range(100):
        A = rng.random((n, n))
        b = rng.random(n)
        A_b[N.index(n)].append((A, b))

# A= A_b[N.index(n)][i][0] b= A_b[N.index(n)][i][1]となる