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
def QR_method(A): #input matrix A
    n= len(A)
    k= 0
    e=1e-6

    while k<10**4:
        k+=1

        if np.max(np.abs(np.tril(A, -1))) < e: break

        Q, R= QR_decomposition(A)
        A = R @ Q

    return np.diag(A).copy(), k


lambda_tildes= [[] for _ in range(len(N))]
W = [[] for _ in range(len(N))]
V = [[] for _ in range(len(N))]
times = [[] for _ in range(len(N))]
iterations = [[] for _ in range(len(N))]

QR_method(A[0][0]) # 最初の呼び出しでNumbaのコンパイルを済ませておく（計算時間に含めないため）

for i in range(len(N)):
    for j in range(100):
        start = time.perf_counter()
        lt, k = QR_method(A[i][j])
        end = time.perf_counter()
        lambda_tildes[i].append(lt)
        iterations[i].append(k)
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


# n, 残差ノルムの最大値, 相対誤差の最大値, 計算時間, 反復回数 を1行ずつ保存（平均しない）
data = []
for i in range(len(N)):
    for j in range(100):
        data.append([N[i], max_residual_norms[i][j], max_rel_errors[i][j], times[i][j], iterations[i][j]])
np.savetxt("result.dat", data, header="n max_residual max_rel_error time iteration", fmt="%.6e")
