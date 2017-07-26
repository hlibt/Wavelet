reset
set term gif animate
set output "animate.gif"
i=0
n=1000
set xrange[0:1]
set yrange[-1:3]
load "animate.gnuplot"
set output
