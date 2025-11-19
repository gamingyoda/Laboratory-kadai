set size ratio -1
set view map
set contour base
unset surface
set xlabel 'x'
set ylabel 'y'
set object 1 rect from 30,0 to 60,40 fc rgb 'gray' fs solid 0.4 border lc rgb 'black'
set cntrparam levels discrete 0,1,2,3,4,5,6,7,8,9,10,20,30
splot 'psi.dat' using 1:2:3 with lines
