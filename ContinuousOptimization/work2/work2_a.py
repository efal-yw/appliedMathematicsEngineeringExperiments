#最急降下法により停留点を求める

def f(x):
    return 1/3 * x**3 - x**2 - 3*x + 5/3

x_k = 1/2
eps = 1e-6

def d(f,x):
    h = 1e-6
    return (f(x+h)-f(x))/h

for k in range(10**5):
    if abs(d(f,x_k)) < eps: break
    d_k = -d(f,x_k)
    t_k = 1/ (k+1)
    x_k += d_k * t_k

print("stationary point:", x_k)