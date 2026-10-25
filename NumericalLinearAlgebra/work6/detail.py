# 課題6 の詳細解析用（QR法は再実行しない）
# result.dat（work6.py の出力）に、各行列の固有値から計算した収束率 (2.16) を付け加えて detail.dat に保存する
import numpy as np
from generateA_b import A, N

A = [[(a + a.T) / 2 for a in row] for row in A]   # work6.py と同じ対称化
d = np.loadtxt("result.dat")                       # n, max_residual, max_rel_error, time, iteration

data = []
for i in range(len(N)):
    for j in range(100):
        w = np.linalg.eigvalsh(A[i][j])
        w = w[np.argsort(-np.abs(w))]               # |λ1| > |λ2| > ... の順に並べる
        ratios = np.abs(w[1:]) / np.abs(w[:-1])     # |λ_{i+1}| / |λ_i|
        m = np.argmax(ratios)
        rho = ratios[m]                             # 収束率 (2.16)
        opposite = 1 if w[m] * w[m+1] < 0 else 0    # 収束率を与えるペアが異符号なら1
        row = d[i*100 + j]
        data.append([N[i], rho, 1 - rho, opposite, row[1], row[2], row[3], row[4], row[3] / row[4]])

np.savetxt("detail.dat", data,
           header="n rho 1-rho opposite_sign max_residual max_rel_error time iteration time_per_iteration",
           fmt="%.6e")
