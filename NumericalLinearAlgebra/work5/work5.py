import numpy as np

rng = np.random.default_rng(seed=0)
n = 50
A = rng.random((n, n))
e = 10**(-6)

# べき乗法を1回だけ行い、λ1の近似値をsigmaとする
x = np.ones(n)
y = A @ x
i = np.argmax(abs(x))
sigma = y[i] / x[i]

# 逆反復法 (A - sigma*I) y = x を解いてλ1に収束させる
B = A - sigma * np.eye(n)
x = y / np.linalg.norm(y) #べき乗法1回目で得たベクトルを初期ベクトルにする
prev_lambda = None #収束判定のためのλ[k-1]
k = 1 #カウンター
while k < 10**4:
    y = np.linalg.solve(B, x)
    i = np.argmax(abs(x))
    my = y[i] / x[i] #(A - sigma*I)^(-1) の固有値の近似
    lambda_1 = sigma + 1 / my
    x = y / np.linalg.norm(y)
    if (prev_lambda is not None and abs(lambda_1 - prev_lambda) / abs(lambda_1) < e): break
    prev_lambda = lambda_1
    k += 1


