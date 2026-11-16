#fxのゼロ点をニュートン法により求める
eps = 1e-6
def f(x):
    return x**3 + 2*x**2 - 5*x - 6

def d(f, x):
    h = 1e-6
    return (f(x+h) - f(x)) / h


def newton(f, x0):
    if abs(f(x0)) < eps: return x0
    x1 = x0 - f(x0)/d(f,x0) 
    return newton(f, x1)

x0s = [-2.5, 0, 2.5] #根拠はmatlibの実際のグラフ

roots = []
for x0 in x0s:
    root = newton(f, x0)
    roots.append(root)

print("Roots:", roots)