set title '(1) Spiral attractor'
set xlabel 'z'
set ylabel 'x'
set xrange[-6:6]
set yrange[-3:3]
plot 'spiral.dat' using 1:2 with lines notitle
