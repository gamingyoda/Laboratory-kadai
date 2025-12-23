set grid
set pm3d map
set title 'Sod 2D: Velocity u (t=0.200, NX=400, NY=400)'
set xlabel 'x'
set ylabel 'y'
set cblabel 'Velocity u'
splot 'sod2d.dat' using 1:2:4 with pm3d notitle
