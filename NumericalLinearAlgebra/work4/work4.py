from generateA_b import A, N
import numpy as np
import time 

for n in range(N):
    
    x = np.ones(n)
    e = 10**(-6)
    prev_my = None #収束判定のためのmy[k-1]

    ans = []*len(N) #結果を収納する
    k=1 #カウンター
    start = time.start()

    while True:
        y = A @ x
        i = np.argmax(abs(x))
        my = y[i][i] / x[i]
        x = y[i][i] / abs(y[i][i])
        end = time.end()
        if (prev_my is not None & (my - prev_my) / my < e): break
        prev_my = my
        k+=1

    #結果の処理
    time_len = end - start
    lambda_1 = my 
    x_1 = x
    real_labmda, real_x = np.linalg.eig(A)
    rel_error = abs(real_x - x_1) / abs(real_x)
    residual_norm = abs(lambda_1 * x_1 - A @ x_1)
    ans[N.index(n)].append([lambda_1, x_1, real_labmda, real_x, rel_error, residual_norm, k, time_len])






    