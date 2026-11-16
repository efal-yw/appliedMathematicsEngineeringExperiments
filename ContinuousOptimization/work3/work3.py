import numpy as np

#(a)
def f(x):
    return x[0]**2 + np.exp(x[0]) + x[1]**4 + x[1]**2 - 2*x[0]*x[1] +3

#(b)
def divf(x):
    h =1e-6
    delx0 = (f([x[0]+h,x[1]])-f([x[0],x[1]]))/h
    delx1 = (f([x[0],x[1]+h])-f([x[0],x[1]]))/h
    return np.array([delx0,delx1])

def Hessef(x):
    h =1e-6
    [delx0delx0, delx1delx0] = (divf(x+np.array([h,0]))-divf(x))/h
    [delx0delx1, delx1delx1] = (divf(x+np.array([0,h]))-divf(x))/h
    H =np.array([[delx0delx0, delx1delx0], [delx0delx1, delx1delx1]])
    return H
