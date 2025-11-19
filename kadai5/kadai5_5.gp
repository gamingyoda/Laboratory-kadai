set xlabel 'x, xi'
set ylabel 'iy, ieta'
set grid
set size ratio -1
set xrange[-2.5:2.5]
set yrange[-2.5:2.5]
set title 'a=1, xi0=-0.2, eta0=0'
set label 1 sprintf('x_c/a = -0.200', -0.2) at 0.6, 2.1
set label 2 sprintf('y_c/a = 0.000', 0) at 0.6, 1.8
set label 3 sprintf('r/a   = 1.200', 1.2) at 0.6, 1.5
plot \
  'circle_5.dat' using 3:4 with lines lw 2 lc rgb 'black' notitle, \
  'joukowski_5.dat' using 3:4 with lines lw 2 lc rgb 'blue'  notitle
