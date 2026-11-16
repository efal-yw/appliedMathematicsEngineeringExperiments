import numpy as np

def x(t) : return t - np.floor(t)

T = 1

def a(k):
    M=10**5
    i = np.arange(-M, M)
    return np.sum(x((T/2) * i/M) * np.cos((2 * np.pi * k / T) * (T/2) * i/M) * (T/2)/M) * 2 / T

def b(k):
    M=10**5
    i = np.arange(-M, M)
    return np.sum(x((T/2) * i/M) * np.sin((2 * np.pi * k / T) * (T/2) * i/M) * (T/2)/M) * 2 / T

def xM(t, M):
    sum = 0
    for i in range(M):
        sum += a(i+1) * np.cos((2 * np.pi * (i+1) / T) * t) + b(i+1) * np.sin((2 * np.pi * (i+1) / T) * t)
    return a(0)/2 + sum


Ms = [10, 10**3, 10**5]

import matplotlib.pyplot as plt

t = np.linspace(-1, 2, 1000)

plt.plot(t, x(t), label="x(t)")
for M in Ms:
    plt.plot(t, xM(t, M), label="M=" + str(M))

plt.xlabel("t")
plt.legend()
plt.show()
