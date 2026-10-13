import numpy as np
import sys
sys.path.append("../..")
from NumericalLinearAlgebra.generateA_b import A_b, N
import time


Ab = A_b.copy()


#LU分解する関数を定義。返り値は (pi, L, U) で、P_pi A = L U となる
#pi は置換を配列で表したもの（P_pi A = A[pi]）
def lu_decomposition(A): #行列Aを引数にとる
    n = A.shape[0]
    A = A.copy()         
    L = np.zeros((n, n))
    pi = np.arange(n)        #最初 pi は恒等置換
    for k in range(n-1):
        #(1) 第k列の下三角成分 a_kk,...,a_nk から絶対値最大のものをピボットに選ぶ
        i = k + np.argmax(np.abs(A[k:, k]))
        #(2) i != k なら A, L の第k行と第i行を入れ替え、pi に左から互換(k i)を掛ける
        if i != k:
            A[[k, i], :] = A[[i, k], :]
            L[[k, i], :] = L[[i, k], :]
            pi[[k, i]] = pi[[i, k]]
        #(4) L の第k列を (0_{k-1}, 1, a_{k+1,k}/a_kk, ..., a_nk/a_kk) とする（(3)で書き換わる前に計算）
        l = A[k+1:, k] / A[k, k]
        L[k, k] = 1
        L[k+1:, k] = l
        #(3) k+1行目以降の各行iから (a_ik / a_kk) * (k行目) を引く
        A[k+1:, :] -= np.outer(l, A[k, :])
    L[n-1, n-1] = 1
    U = A                    #最後の A を U とする
    return pi, L, U


#LU分解を用いて Ax = b を解く関数を定義。返り値は解ベクトルx
def lu_solve(A, b): #行列A、ベクトルbを引数にとる
    n = len(b)
    #(1) P_pi A = L U を満たす L, U, pi を求める
    pi, L, U = lu_decomposition(A)
    #(2) Ly = P_pi b を前進代入で解く：y_i = b_pi(i) - Σ_{j=1}^{i-1} l_ij y_j
    y = np.zeros(n)
    for i in range(n):
        y[i] = b[pi[i]] - L[i, :i] @ y[:i]
    #(3) Ux = y を後退代入で解く：x_i = (y_i - Σ_{j=i+1}^{n} u_ij x_j) / u_ii
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - U[i, i+1:] @ x[i+1:]) / U[i, i]
    return x


X_n =[] #解いたやつ
X_n_performance_time=[]
X_n_res = [] #残差ノルム　X_n_res[n.index]　各nで平均をとっている
X_n_rel_res = [] #相対誤差
for n_index in range(len(Ab)):
    X_n_list = []
    X_n_performance_time_i = []
    X_n_res_i = []
    X_n_rel_res_i = []
    for i in range(100):
        A = Ab[n_index][i][0]
        b = Ab[n_index][i][1]
        start = time.perf_counter()
        X_n_i = lu_solve(A, b)
        end = time.perf_counter()
        X_n_list.append(X_n_i)
        X_n_performance_time_i.append(end - start)

        r = b - A @ X_n_i                         # 残差ベクトル
        res_norm = np.linalg.norm(r)           # 2ノルム
        x_true = np.linalg.solve(A, b)         # 真の解の代わりにnumpyの解を使う
        rel_err = np.linalg.norm(X_n_i - x_true) / np.linalg.norm(x_true)  # 相対誤差

        X_n_res_i.append(res_norm)
        X_n_rel_res_i.append(rel_err)
    X_n_performance_time.append(np.average(X_n_performance_time_i))
    X_n_res.append(np.average(X_n_res_i))
    X_n_rel_res.append(np.average(X_n_rel_res_i))
    X_n.append(X_n_list)

# n, 残差ノルム, 相対誤差, 計算時間 を列にして保存
data = np.column_stack((N, X_n_res, X_n_rel_res, X_n_performance_time))
np.savetxt("result.dat", data, header="n residual relative_error time", fmt="%.6e")

