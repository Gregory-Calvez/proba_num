set nokey
set xlabel "Inverse of number of steps"
plot "plot.dat" using 1:2:3 with yerrorlines, \
"plot.dat" using 1:4:5 with yerrorlines, \
6.144with lines lt 3