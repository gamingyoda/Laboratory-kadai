set xlabel 'x, xi'
set ylabel 'iy, ieta'
set grid
set size ratio -1
set xrange[-2.5:2.5]
set yrange[-2.5:2.5]
set title 'a=1, xi0=0, eta0=0.3'
set label 1 sprintf('x_c/a = 0.000', 0) at 0.6, 2.1
set label 2 sprintf('y_c/a = 0.300', 0.3) at 0.6, 1.8
set label 3 sprintf('r/a   = 1.044', 1.04403) at 0.6, 1.5
plot \
  'circle_3.dat' using 3:4 with lines lw 2 lc rgb 'black' notitle, \
  'joukowski_3.dat' using 3:4 with lines lw 2 lc rgb 'blue'  notitle
