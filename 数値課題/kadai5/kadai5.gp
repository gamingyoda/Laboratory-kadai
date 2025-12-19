set xlabel 'x'
set ylabel 'y'
set grid
set key outside
set size ratio -1
plot \
  'joukowski_1.dat' using 1:2 with lines lw 2 lc rgb 'black' title 'case 1', \
  'joukowski_2.dat' using 1:2 with lines lw 2 lc rgb 'red'   title 'case 2', \
  'joukowski_3.dat' using 1:2 with lines lw 2 lc rgb 'blue'  title 'case 3', \
  'joukowski_4.dat' using 1:2 with lines lw 2 lc rgb 'green' title 'case 4', \
  'joukowski_5.dat' using 1:2 with lines lw 2 lc rgb 'magenta' title 'case 5', \
  'joukowski_6.dat' using 1:2 with lines lw 2 lc rgb 'cyan' title 'case 6', \
  'joukowski_7.dat' using 1:2 with lines lw 2 lc rgb 'orange' title 'case 7', \
  'joukowski_8.dat' using 1:2 with lines lw 2 lc rgb 'brown' title 'case 8', \
  'joukowski_9.dat' using 1:2 with lines lw 2 lc rgb 'violet' title 'case 9'
