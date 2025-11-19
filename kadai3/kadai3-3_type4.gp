set title '(4) Sparrow attractor'
set xlabel 'y'
set ylabel 'x'
set xrange[-2:2]
set yrange[-2:2]
plot 'sparrow.dat' using 1:2 with lines notitle
