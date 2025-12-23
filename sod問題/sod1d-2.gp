set grid
set multiplot layout 3,1 title 'Sod 1D: Roe-FDS + TVD RK2 (t=0.200, NX=1000)'
set xlabel 'x'
set ylabel 'rho'
plot 'sod1d.dat' using 1:2 with lines title 'rho'
set xlabel 'x'
set ylabel 'u'
plot 'sod1d.dat' using 1:3 with lines title 'u'
set xlabel 'x'
set ylabel 'p'
plot 'sod1d.dat' using 1:4 with lines title 'p'
unset multiplot
