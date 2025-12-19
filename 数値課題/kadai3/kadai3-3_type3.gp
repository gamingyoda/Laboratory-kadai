set title '(3) Double screw attractor'
set xlabel 'y'
set ylabel 'z'
set xrange[-70:70]
set yrange[-70:70]
plot 'dscrew.dat' using 1:2 with lines notitle
