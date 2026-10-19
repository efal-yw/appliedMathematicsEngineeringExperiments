import numpy as np
from generateA_b import A, N
n=50 #THIS IS A lIE


def QR_method(A): #input matrix A
    n= len(A)
    k= 0
    e=10e-6


    def QR_decomposition(X):
        n= len(X)
        R= np.zeros((n,n))
        Q= np.zeros((n,n))

        for  j in range(n):
            v = X[:,j].astype(float)
            for i in range(j):
                R[i,j] =  Q[:,i] @ X[:,j]
                v = v-R[i,j]*Q[:,i]
            R[j,j] = np.linalg.norm(v)
            Q[:,j] = v/R[j,j]
        return Q, R

    
    while k<10**4:
        k+=1

        convergence_check_counter= 0
        for i in range(n):
            for j in range(i):
                if abs(A[i][j])>e: convergence_check_counter+=1
        if convergence_check_counter==0: break 

        Q, R= QR_decomposition(A)
        A = R @ Q

    return np.array([A[i,i] for i in range(n)])


results = []
for i in range(len(N)):
    for j in range(100):
        results.append(QR_method(A[i][j]))