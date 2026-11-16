import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return [x_**3 + 2*x_**2 -5*x_ -6 for x_ in x]

# 課題１
x = np.arange(-10,10,0.1)
y = f(x)
plt.plot(x, y, "red")
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.axhline(0, c="black")
plt.axvline(0, c="black")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Graph of y=f(x)")
plt.show()

def dichotomy(f,a,b):
    eps = 1e-6

    

dichotomy(f, 1, 3)