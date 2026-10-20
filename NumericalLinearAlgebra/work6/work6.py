import numpy as np
from generateA_b import A, N
import time

A = [[(a+a.T)/2 for a in row] for row in A]
def QR_method(A): #input matrix A
    n= len(A)
    k= 0
    e=1e-6


    def QR_decomposition(X):
        n= len(X)
        R= np.zeros((n,n))
        Q= np.zeros((n,n))

        for  j in range(n):
            v = X[:,j].astype(float)
            for i in range(j):
                R[i,j] =  Q[:,i] @ v
                v = v-R[i,j]*Q[:,i]
            R[j,j] = np.linalg.norm(v)
            Q[:,j] = v/R[j,j]
        return Q, R

    
    while k<10**4:
        k+=1

        if np.max(np.abs(np.tril(A, -1))) < e: break

        Q, R= QR_decomposition(A)
        A = R @ Q

    return np.array([A[i,i] for i in range(n)])


lambda_tildes= [[] for _ in range(len(N))]
W = [[] for _ in range(len(N))]
V = [[] for _ in range(len(N))]
times = [[] for _ in range(len(N))]

for i in range(len(N)):
    for j in range(100):
        start = time.perf_counter()
        lambda_tildes[i].append(QR_method(A[i][j]))
        end = time.perf_counter()
        il_tempo = end - start
        times[i].append(il_tempo)
        lambda_tildes[i][j] = lambda_tildes[i][j][np.argsort(-np.abs(lambda_tildes[i][j]))]
        w, v = np.linalg.eigh(A[i][j])
        ids = np.argsort(-np.abs(w))
        W[i].append(w[ids])
        V[i].append(v[:, ids])


max_residual_norms = [[] for _ in range(len(N))]
max_rel_errors = [[] for _ in range(len(N))]
for i in range(len(N)):
    for j in range(100):
        a, lt, w, v = A[i][j], lambda_tildes[i][j], W[i][j], V[i][j]
        max_residual_norms[i].append(max(np.linalg.norm((a - lt[k]*np.eye(N[i])) @ v[:, k]) for k in range(N[i])))
        max_rel_errors[i].append(np.max(np.abs(w - lt) / np.abs(w)))
