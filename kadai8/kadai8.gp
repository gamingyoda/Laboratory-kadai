set size ratio -1
set view map
set contour base
unset surface
set cntrparam levels discrete 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0
set object 1 rect from 0,0 to 4,1 fc rgb 'gray' fs solid 0.4 border lc rgb 'black'
splot 'field.dat' using 1:2:3 with lines
