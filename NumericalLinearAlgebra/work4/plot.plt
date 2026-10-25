set terminal pngcairo size 800,600 font "Arial,12"
set logscale xy                 # 両対数
set xlabel "n"
set grid
set key top left
set format y "10^{%L}"          # 縦軸を 10^x の形で表示

# 平均の線 (smooth unique で同じ n の値を平均する)
avg = 'with lines lc rgb "red" lw 2 title "average"'

# 残差ノルム (2.12)
set output "residual.png"
set ylabel "residual norm ||A v~_1 - {/Symbol l}~_1 v~_1||_2"
plot "result.dat" using 1:2 with points pt 7 ps 0.5 title "residual norm", \
     "" using 1:2 smooth unique @avg

# λ~1 の相対誤差 (2.10)
set output "relative_error_lambda.png"
set ylabel "relative error of {/Symbol l}~_1"
plot "result.dat" using 1:3 with points pt 7 ps 0.5 title "relative error ({/Symbol l}_1)", \
     "" using 1:3 smooth unique @avg

# v~1 の相対誤差 (2.11)
set output "relative_error_x.png"
set ylabel "relative error of v~_1"
plot "result.dat" using 1:4 with points pt 7 ps 0.5 title "relative error (v_1)", \
     "" using 1:4 smooth unique @avg

# 計算時間 (n^2 の参考線つき)
set output "time.png"
set ylabel "time / s"
stats "result.dat" using 5 every ::0::99 nooutput   # n=50 の100個分の time
t0 = STATS_mean
plot "result.dat" using 1:5 with points pt 7 ps 0.5 title "time", \
     "" using 1:5 smooth unique @avg, \
     t0*(x/50.0)**2 with lines dt 2 title "O(n^2)"

# 反復回数
set output "iteration.png"
unset logscale y
set format y "%g"
set ylabel "iterations"
plot "result.dat" using 1:6 with points pt 7 ps 0.5 title "iterations", \
     "" using 1:6 smooth unique @avg

unset output
