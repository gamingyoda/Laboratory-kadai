set xlabel 'x/a'
set ylabel 'iy/a'
set grid
set xrange [-5:5]
set yrange [-5:5]
set size ratio -1
set label 1 sprintf('{/Symbol x}/a = -0.100, {/Symbol h}/a = 0.250, r/a = 1.128', -0.1, 0.25, 1.12805141726785) at graph 0.98,0.10 right
set label 2 sprintf('{/Symbol G}/(2{/Symbol p} a U_0) = -1.052, {/Symbol a} = 15.0 deg', -1.05236481237008, 15) at graph 0.98,0.05 right
plot \
  'airfoil.dat' using ($1/1):($2/1) with lines lt -1 notitle, \
  'stream_00.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_01.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_02.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_03.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_04.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_05.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_06.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_07.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_08.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_09.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_10.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_11.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_12.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_13.dat' using ($1/1):($2/1) with lines lw 1 notitle, \
  'stream_14.dat' using ($1/1):($2/1) with lines lw 1 notitle
