set size ratio -1
set view map
set contour base
unset surface
set cntrparam levels incremental 0,0.2,4.0
set object 1 rect from 0,0 to 4,1 fc rgb 'gray' fs solid 0.4 border lc rgb 'black'
splot 'field.dat' using 1:2:3 with lines notitle
