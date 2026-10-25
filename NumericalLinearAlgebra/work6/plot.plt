set terminal pngcairo size 800,600 font "Arial,12"
set logscale xy                 # 両対数
set xlabel "n"
set grid
set key top left
set format y "10^{%L}"          # 縦軸を 10^x の形で表示

# 平均の線 (smooth unique で同じ n の値を平均する)
avg = 'with lines lc rgb "red" lw 2 title "average"'

# 固有方程式の残差ノルムの最大値
set output "residual.png"
set ylabel "max_i ||(A - {/Symbol l}~_i I) v_i||_2"
plot "result.dat" using 1:2 with points pt 7 ps 0.5 title "max residual norm", \
     "" using 1:2 smooth unique @avg

# 相対誤差の最大値
set output "relative_error.png"
set ylabel "max_i |{/Symbol l}_i - {/Symbol l}~_i| / |{/Symbol l}_i|"
plot "result.dat" using 1:3 with points pt 7 ps 0.5 title "max relative error", \
     "" using 1:3 smooth unique @avg

# 計算時間 (n^3 の参考線つき)
set output "time.png"
set ylabel "time / s"
stats "result.dat" using 4 every ::0::99 nooutput   # n=10 の100個分の time
t0 = STATS_mean
plot "result.dat" using 1:4 with points pt 7 ps 0.5 title "time", \
     "" using 1:4 smooth unique @avg, \
     t0*(x/10.0)**3 with lines dt 2 title "O(n^3)"

# 反復回数
set output "iteration.png"
set ylabel "iterations"
plot "result.dat" using 1:5 with points pt 7 ps 0.5 title "iterations", \
     "" using 1:5 smooth unique @avg

unset output
