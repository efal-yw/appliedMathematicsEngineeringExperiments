#fxのゼロ点を二分法により求める

def f(x):
    return x**3 + 2*x**2 - 5*x - 6


def dichotomy(f,a,b): #a<b and f(a)*f(b)<0
    def f_(x):
        return f(x) * f(b) // abs(f(b))
    eps = 1e-6
    c= (a+b)/2

    while abs(f(c)) > eps:
        if f_(c) < 0:
            a = c
        else:
            b = c
        c = (a+b)/2
    return c


a_bs = [(-5, -2.5), (-2.5, -0), (1, 2.5)] #根拠はmatlibの実際のグラフ

roots = []
for a, b in a_bs:
    root = dichotomy(f, a, b)
    roots.append(root)

print("Roots:", roots)