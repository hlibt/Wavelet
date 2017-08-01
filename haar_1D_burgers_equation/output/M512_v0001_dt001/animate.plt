reset
set term gif animate
set output "animate.gif"
i=0
n=1000
set xrange[0:1.25]
set yrange[0:1]
load "animate.gnuplot"
set output
