set xlabel 'dimensionless distance x'
set ylabel 'dimensionless temperature u(x,t)'
set grid
plot \
  'heat_t0.dat' using 1:2 with lines lw 2 lc rgb 'black' title 't=0', \
  'heat_t1000.dat' using 1:2 with lines lw 1 lc rgb 'blue' title 't=0.04', \
  'heat_t2000.dat' using 1:2 with lines lw 1 lc rgb 'blue' title 't=0.08', \
  'heat_t3000.dat' using 1:2 with lines lw 1 lc rgb 'blue' title 't=0.12', \
  'heat_t4000.dat' using 1:2 with lines lw 1 lc rgb 'blue' title 't=0.16', \
  'heat_t5000.dat' using 1:2 with lines lw 1 lc rgb 'blue' title 't=0.2', \
  'heat_tend.dat' using 1:2 with lines lw 2 lc rgb 'red' title 't=end'
