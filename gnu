set nokey
set xlabel "Inverse of number of steps"
set xrange [0.005:0.205]
plot "plot.dat" using 1:2:3 with yerrorlines, \
"plot.dat" using 1:4:5 with yerrorlines, \
4.000000 with lines lt 3