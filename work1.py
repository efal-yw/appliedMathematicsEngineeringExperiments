import numpy as np
from generateA_b import A_b
import time


Ab = A_b.copy()


#前進消去する関数を定義。返り値は拡大係数行列A_b
def forward_eraser(A, b): #行列A、ベクトルbを引数にとる
    n = len(b)
    A_b = np.column_stack((A, b)) #Aとbを結合して拡大係数行列を作る
    for k in range(n-1):
        A_kMaxIndex = k + np.argmax(np.abs(A_b[k:, k]))
        if A_kMaxIndex != k:
            A_b[A_kMaxIndex, :], A_b[k, :] = A_b[k, :].copy(), A_b[A_kMaxIndex, :].copy()  
        for i in range(k+1, n):
            A_b[i, :] = A_b[i, :] - (A_b[i, k] / A_b[k, k]) * A_b[k, :]
    return A_b

# A= A_b[N.index(n)][i][0] b= A_b[N.index(n)][i][1]となる

#漸化式（２）

def gaussian_elimination_i(i, x, Ab_prime):   #i番目の解を求める（x[i+1:]は求め済み）
    n = Ab_prime.shape[0]
    sum_aijxj = 0
    for j in range(i+1, n):
        sum_aijxj += Ab_prime[i][j] * x[j]
    x_i = (1/Ab_prime[i][i]) * (Ab_prime[i][n] - sum_aijxj)
    return x_i


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
        start = time.perf_counter()
        A = Ab[n_index][i][0]
        b = Ab[n_index][i][1]
        Ab_prime = forward_eraser(A, b)
        #Ab_primeを準備した。あとは漸化式に入れるだけ
        real_n = Ab_prime.shape[0]
        X_n_i = np.zeros(real_n)
        for k in range(real_n):
            row = real_n - 1 - k   #下の行から順に解く
            X_n_i[row] = gaussian_elimination_i(row, X_n_i, Ab_prime)
        X_n_list.append(X_n_i)
        end= time.perf_counter()
        performance_time_i = end - start
        X_n_performance_time_i.append(performance_time_i)

        r = b - A @ X_n_i                         # 残差ベクトル（@ は行列×ベクトルの積）
        res_norm = np.linalg.norm(r)           # 2ノルム
        x_true = np.linalg.solve(A, b)         # 真の解の代わりにnumpyの解を使う
        rel_err = np.linalg.norm(X_n_i - x_true) / np.linalg.norm(x_true)  # 相対誤差

        X_n_res_i.append(res_norm)
        X_n_rel_res_i.append(rel_err)
    X_n_performance_time.append(np.average(X_n_performance_time_i))
    X_n_res.append(np.average(X_n_res_i))
    X_n_rel_res.append(np.average(X_n_rel_res_i))
    X_n.append(X_n_list)

real_X_n = []  #numpyにより解いた解
for n_index in range(len(Ab)):
    real_X_n_list = []
    for i in range(100):
        A = Ab[n_index][i][0]
        b = Ab[n_index][i][1]
        real_X_n_list.append(np.linalg.solve(A, b))
    real_X_n.append(real_X_n_list)


