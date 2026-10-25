set terminal pngcairo size 800,600 font "Arial,12"
set logscale xy                 # 両対数
set xlabel "n"
set grid
set key top left
set format y "10^{%L}"          # 縦軸を 10^x の形で表示

# 残差ノルム (2.5)
set output "residual.png"
set ylabel "residual norm ||b - A x~||_2"
plot "result.dat" using 1:2 with linespoints pt 7 title "residual norm"

# 相対誤差 (2.4)
set output "relative_error.png"
set ylabel "relative error"
plot "result.dat" using 1:3 with linespoints pt 7 title "relative error"

# 計算時間 (n^3 の参考線つき)
set output "time.png"
set ylabel "time / s"
stats "result.dat" using 4 every ::0::0 nooutput   # 最初のデータ行 (n=100) の time
t0 = STATS_min
plot "result.dat" using 1:4 with linespoints pt 7 title "time", \
     t0*(x/100.0)**3 with lines dt 2 title "O(n^3)"

unset output
