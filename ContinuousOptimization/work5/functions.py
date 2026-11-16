import numpy as np

y=[1.5, 2.25, 2.625]
def fis(x):
    fis = [y[i] - x[0] * (1- x[1]**(i+1)) for i in range(3)]
    return fis

def f(x):
    sum =0
    for fi in fis(x):
        sum+=fi**2
    return sum


def divf(x):
    h =1e-6
    delx0 = (f([x[0]+h,x[1]])-f([x[0]-h,x[1]]))/(2*h)
    delx1 = (f([x[0],x[1]+h])-f([x[0],x[1]-h]))/(2*h)
    return np.array([delx0,delx1])

def Hessef(x):
    h =1e-6
    [delx0delx0, delx1delx0] = (divf(x+np.array([h,0]))-divf(x-np.array([h,0])))/(2*h)
    [delx0delx1, delx1delx1] = (divf(x+np.array([0,h]))-divf(x-np.array([0,h])))/(2*h)
    H =np.array([[delx0delx0, delx1delx0], [delx0delx1, delx1delx1]])
    return H
