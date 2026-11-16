#ニュートン法により停留点を求める

def f(x):
    return 1/3 * x**3 - x**2 - 3*x + 5/3

x_k = 5
eps = 1e-6
t_k = 1

def d(f,x):
    h = 1e-6
    return (f(x+h)-f(x))/h

def d2(f, x):
    h = 1e-6
    return (d(f,x+h)-d(f,x))/h

for k in range(10**5):
    if abs(d(f,x_k)) < eps: break
    d_k = -d(f,x_k)/d2(f,x_k)
    x_k += d_k * t_k

print("stationary point:", x_k)