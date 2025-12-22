set grid
set pm3d map
set title 'Sod 2D: Roe-FDS + TVD RK2 (t=0.200, NX=1000, NY=1000)'
set xlabel 'x'
set ylabel 'y'
set cblabel 'Density rho'
splot 'sod2d.dat' using 1:2:3 with pm3d notitle
