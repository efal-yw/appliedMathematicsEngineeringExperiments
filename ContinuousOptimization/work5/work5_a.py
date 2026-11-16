from functions import f, divf, Hessef
import numpy as np

x_0 = np.array([2.0,0.0])
xi = 1e-4
rho = 0.5
tau_0 = 1
eps = 1e-5

def tau(x_k,d_k):
    tau = tau_0
    while f(x_k + tau*d_k) > f(x_k) + xi * tau * (d_k @ divf(x_k)):
        tau = rho * tau
    return tau


def steepest_descent_method():
    x_k = x_0.copy()
    for k in range(10**8):
        if np.linalg.norm(divf(x_k)) <= eps:
            return k, x_k
        d_k = -divf(x_k)
        x_k += tau(x_k, d_k) * d_k

k, x_k = steepest_descent_method()

print(k, x_k)