# 課題3: 消去法 (work1) と LU分解による解法 (work2) の比較
set terminal pngcairo size 800,600 font "Arial,12"
set logscale xy
set xlabel "n"
set grid
set key top left
set format y "10^{%L}"

E = "../work1/result.dat"
L = "../work2/result.dat"

set output "residual_compare.png"
set ylabel "residual norm ||b - A x~||_2"
plot E using 1:2 with linespoints pt 7 title "Gaussian elimination", \
     L using 1:2 with linespoints pt 5 dt 2 title "LU decomposition"

set output "relative_error_compare.png"
set ylabel "relative error"
plot E using 1:3 with linespoints pt 7 title "Gaussian elimination", \
     L using 1:3 with linespoints pt 5 dt 2 title "LU decomposition"

set output "time_compare.png"
set ylabel "time [s]"
stats E using 4 every ::0::0 nooutput
t0 = STATS_min
plot E using 1:4 with linespoints pt 7 title "Gaussian elimination", \
     L using 1:4 with linespoints pt 5 dt 2 title "LU decomposition", \
     t0*(x/100.0)**3 with lines dt 3 title "O(n^3)"

unset output
