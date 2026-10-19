from generateA_b import A, N
import numpy as np
import time


#LU分解する関数を定義（work2と同じ）。返り値は (pi, L, U) で、P_pi A = L U となる
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


ans = [[] for _ in N] #結果を収納する
e = 10**(-6)
for n in N:
    for _ in range(100):
        A_ = A[N.index(n)][_]
        start = time.perf_counter()

        # べき乗法を1回だけ行い、λ1の近似値をsigmaとする
        x = np.ones(n)
        y = A_ @ x
        i = np.argmax(abs(x))
        sigma = y[i] / x[i]

        # 逆反復法 (A - sigma*I) y = x を解いてλ1に収束させる
        B = A_ - sigma * np.eye(n)
        pi, L, U = lu_decomposition(B) #初めに1回だけLU分解しておく
        x = y / np.linalg.norm(y) #べき乗法1回目で得たベクトルを初期ベクトルにする
        prev_lambda = None #収束判定のためのλ[k-1]
        k = 1 #カウンター
        while k < 10**4:
            #Lz = P_pi x を前進代入で解く
            z = np.zeros(n)
            for j in range(n):
                z[j] = x[pi[j]] - L[j, :j] @ z[:j]
            #Uy = z を後退代入で解く
            y = np.zeros(n)
            for j in range(n-1, -1, -1):
                y[j] = (z[j] - U[j, j+1:] @ y[j+1:]) / U[j, j]
            i = np.argmax(abs(x))
            my = y[i] / x[i] #(A - sigma*I)^(-1) の固有値の近似
            lambda_1 = sigma + 1 / my
            x = y / np.linalg.norm(y)
            if (prev_lambda is not None and abs(lambda_1 - prev_lambda) / abs(lambda_1) < e): break
            prev_lambda = lambda_1
            k += 1
        end = time.perf_counter()

        #結果の処理
        time_len = end - start
        x_1 = x
        real_lambdas, real_xs = np.linalg.eig(A_)
        real_lambda = real_lambdas[np.argmax(real_lambdas)]
        real_x = real_xs[:, np.argmax(real_lambdas)]
        if (real_x @ x_1).real < 0: real_x = -real_x
        rel_error_x = np.linalg.norm(real_x - x_1) / np.linalg.norm(real_x)
        rel_error_lambda = abs(real_lambda - lambda_1) / abs(real_lambda)
        residual_norm = np.linalg.norm(lambda_1 * x_1 - A_ @ x_1)
        ans[N.index(n)].append([lambda_1, x_1, real_lambda, real_x, rel_error_lambda, rel_error_x, residual_norm, k, time_len])


data = []
for n in N:
    for a in ans[N.index(n)]:
        data.append([n, a[6], a[4], a[5], a[8], a[7]])
np.savetxt("result.dat", data, header="n residual rel_error_lambda rel_error_x time iteration", fmt="%.6e")
