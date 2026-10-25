# 課題6 詳細グラフ（detail.dat を使用, 先に detail.py を実行）
set terminal pngcairo size 800,600 font "Arial,12"
set grid
set key top left
set logscale xy
set format x "10^{%L}"
set format y "10^{%L}"

eps = 1e-6

# 反復回数と収束率 (2.16) の関係
set output "iteration_vs_rho.png"
set key bottom left box opaque
set xlabel "1 - {/Symbol r}   ({/Symbol r} = max_i |{/Symbol l}_{i+1}| / |{/Symbol l}_i|)"
set ylabel "iterations"
set xrange [1e-4:0.3]
set yrange [1:2e4]
plot "detail.dat" using ($1==10?$3:1/0):8 with points pt 7 ps 0.6 title "n = 10", \
     "" using ($1==20?$3:1/0):8 with points pt 5 ps 0.6 title "n = 20", \
     "" using ($1==40?$3:1/0):8 with points pt 9 ps 0.6 title "n = 40", \
     "" using ($1==80?$3:1/0):8 with points pt 11 ps 0.6 title "n = 80", \
     log(eps)/log(1-x) with lines dt 2 lc rgb "black" title "log {/Symbol e} / log {/Symbol r}  ({/Symbol e} = 10^{-6})", \
     1e4 with lines dt 3 lc rgb "gray" title "upper limit 10^4"

# 相対誤差の最大値と反復回数の関係
set output "rel_error_vs_iteration.png"
set key top left box opaque
set xlabel "iterations"
set ylabel "max_i |{/Symbol l}_i - {/Symbol l}~_i| / |{/Symbol l}_i|"
set autoscale x
set autoscale y
set xrange [50:2e4]
plot "detail.dat" using ($1==10?$8:1/0):6 with points pt 7 ps 0.6 title "n = 10", \
     "" using ($1==20?$8:1/0):6 with points pt 5 ps 0.6 title "n = 20", \
     "" using ($1==40?$8:1/0):6 with points pt 9 ps 0.6 title "n = 40", \
     "" using ($1==80?$8:1/0):6 with points pt 11 ps 0.6 title "n = 80"

# 1反復あたりの計算時間
set output "time_per_iteration.png"
set key top left nobox
set xlabel "n"
set ylabel "time / iterations [s]"
set autoscale x
set autoscale y
set xrange [8:100]
set format x "%g"
stats "detail.dat" using 9 every ::0::99 nooutput   # n=10 の100個分
t0 = STATS_mean
avg = 'with lines lc rgb "red" lw 2 title "average"'
plot "detail.dat" using 1:9 with points pt 7 ps 0.5 title "time / iterations", \
     "" using 1:9 smooth unique @avg, \
     t0*(x/10.0)**3 with lines dt 2 title "O(n^3)"

unset output
