from generateA_b import A, N
import numpy as np
import time 

ans = [[] for _ in N] #結果を収納する
e = 10**(-6)
for n in N:
    for _ in range(100):
        A_ = A[N.index(n)][_]
        x = np.ones(n)
        prev_my = None #収束判定のためのmy[k-1]
        k=1 #カウンター
        start = time.perf_counter()

        while k<10**4:
            y = A_ @ x
            i = np.argmax(abs(x))
            my = y[i] / x[i]
            x = y / np.linalg.norm(y)
            if (prev_my is not None and abs(my - prev_my) / abs(my) < e): break
            prev_my = my
            k+=1
        end = time.perf_counter()

        #結果の処理
        time_len = end - start
        lambda_1 = my 
        x_1 = x
        real_lambdas, real_xs = np.linalg.eig(A_)
        real_lambda = real_lambdas[np.argmax(real_lambdas)]
        real_x = real_xs[:, np.argmax(real_lambdas)]
        if (real_x @ x_1).real < 0: real_x = -real_x
        rel_error_x = np.linalg.norm(real_x - x_1) / np.linalg.norm(real_x)
        rel_error_lambda = abs(real_lambda - lambda_1) / abs(real_lambda)
        residual_norm = np.linalg.norm(lambda_1 * x_1 - A_ @ x_1)
        ans[N.index(n)].append([lambda_1, x_1, real_lambda, real_x, rel_error_lambda, rel_error_x, residual_norm, k, time_len])






    
# n, 残差ノルム, λ1の相対誤差, v1の相対誤差, 計算時間, 反復回数 を1行ずつ保存（平均しない）
data = []
for n in N:
    for a in ans[N.index(n)]:
        data.append([n, a[6], a[4], a[5], a[8], a[7]])
np.savetxt("result.dat", data, header="n residual rel_error_lambda rel_error_x time iteration", fmt="%.6e")
