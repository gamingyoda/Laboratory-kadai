set title '(2) Double scroll attractor'
set xlabel 'z'
set ylabel 'x'
set xrange[-6:6]
set yrange[-3:3]
plot 'dscroll.dat' using 1:2 with lines notitle
