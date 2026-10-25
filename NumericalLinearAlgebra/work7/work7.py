import numpy as np
from generateA_b import A, N
import time
from numba import njit

A = [[np.ascontiguousarray((a+a.T)/2) for a in row] for row in A]


@njit
def QR_decomposition(X):
    n= len(X)
    R= np.zeros((n,n))
    Q= np.zeros((n,n))

    for  j in range(n):
        v = X[:,j].astype(np.float64)
        for i in range(j):
            R[i,j] =  np.dot(Q[:,i], v)
            v = v-R[i,j]*Q[:,i]
        R[j,j] = np.linalg.norm(v)
        Q[:,j] = v/R[j,j]
    return Q, R


@njit
def QR_forLambdaNN(A): #input matrix A
    n= len(A)
    k= 0
    e=1e-6
    m=10**4

    while k<m:
        sigma = A[n-1,n-1]
        A =  A - sigma * np.eye(n)
        k+=1
        if np.max(np.abs(A[n-1,:n-1])) < e: 
            A += sigma * np.eye(n)
            break
        Q, R= QR_decomposition(A)
        A = R @ Q
        A += sigma * np.eye(n)

    return A


def QR_speedy(A):
    results = []
    n=len(A)
    for _ in range(n):
        if n==1: 
            results.append(A[0,0]) 
            break
        A= QR_forLambdaNN(A)
        lambda_n = A[-1,-1].copy()
        results.append(lambda_n)
        A = A[:n-1,:n-1]
        n=len(A)
    return results




