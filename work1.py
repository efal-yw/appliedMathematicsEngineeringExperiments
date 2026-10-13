import numpy as np
from generateA_b import A_b

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

def gaussian_elimination_i(i, n_index, Ab_prime):   #i番目の解を求める
    sum_aijxj = 0
    for j in range(i+1,(2**n_index)*100):
      if j>=n: break
      sum_aijxj += Ab_prime[i][j]*gaussian_elimination_i(j, n_index, Ab_prime)
    x_i = (1/Ab_prime[i][i]) * (Ab_prime[(2**n_index)*100] - sum_aijxj)
    return x_i


X_n =[]

for n in range(len([100, 200, 400, 800])):
    for i in range(100):
        A = Ab[n][i][0]
        b = Ab[n][i][1]
        Ab_prime = forward_eraser(A, b)
        #Ab_primeを準備した。あとは漸化式に入れるだけ
        real_n = Ab_prime[(2**n_index)*100]
        _X_n_i=[]
        X_n_i=[]
        for k in range(real_n):
            _X_n_i.append(gaussian_elimination_i(real_n - k))
        for k in range(real_n):
            X_n_i.append(_X_n_i[real_n - k])
    X_n.append(X_n_i)



