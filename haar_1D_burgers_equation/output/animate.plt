reset
set term gif animate
set output "animate.gif"
i=0
n=100
set xrange[0:1]
set yrange[0:1]
load "animate.gnuplot"
set output
