set title '(3) dx/dt = y + x*sqrt(x**2+y**2)*(x**2+y**2-1)**2, dy/dt = -x + y*sqrt(x**2+y**2)*(x**2+y**2-1)**2'
set xrange[-2:2]
set yrange[-2:2]
plot \
'kadai3-2_sys3.dat' using 2:3 with lines linecolor 1 title 'r=0.50', \
'kadai3-2_sys3.dat' using 4:5 with lines linecolor 1 notitle, \
'kadai3-2_sys3.dat' using 6:7 with lines linecolor 1 notitle, \
'kadai3-2_sys3.dat' using 8:9 with lines linecolor 1 notitle, \
'kadai3-2_sys3.dat' using 10:11 with lines linecolor 1 notitle, \
'kadai3-2_sys3.dat' using 12:13 with lines linecolor 1 notitle, \
'kadai3-2_sys3.dat' using 14:15 with lines linecolor 1 notitle, \
'kadai3-2_sys3.dat' using 16:17 with lines linecolor 1 notitle, \
'kadai3-2_sys3.dat' using 18:19 with lines linecolor 2 title 'r=0.90', \
'kadai3-2_sys3.dat' using 20:21 with lines linecolor 2 notitle, \
'kadai3-2_sys3.dat' using 22:23 with lines linecolor 2 notitle, \
'kadai3-2_sys3.dat' using 24:25 with lines linecolor 2 notitle, \
'kadai3-2_sys3.dat' using 26:27 with lines linecolor 2 notitle, \
'kadai3-2_sys3.dat' using 28:29 with lines linecolor 2 notitle, \
'kadai3-2_sys3.dat' using 30:31 with lines linecolor 2 notitle, \
'kadai3-2_sys3.dat' using 32:33 with lines linecolor 2 notitle, \
'kadai3-2_sys3.dat' using 34:35 with lines linecolor 3 title 'r=1.00', \
'kadai3-2_sys3.dat' using 36:37 with lines linecolor 3 notitle, \
'kadai3-2_sys3.dat' using 38:39 with lines linecolor 3 notitle, \
'kadai3-2_sys3.dat' using 40:41 with lines linecolor 3 notitle, \
'kadai3-2_sys3.dat' using 42:43 with lines linecolor 3 notitle, \
'kadai3-2_sys3.dat' using 44:45 with lines linecolor 3 notitle, \
'kadai3-2_sys3.dat' using 46:47 with lines linecolor 3 notitle, \
'kadai3-2_sys3.dat' using 48:49 with lines linecolor 3 notitle, \
'kadai3-2_sys3.dat' using 50:51 with lines linecolor 4 title 'r=1.10', \
'kadai3-2_sys3.dat' using 52:53 with lines linecolor 4 notitle, \
'kadai3-2_sys3.dat' using 54:55 with lines linecolor 4 notitle, \
'kadai3-2_sys3.dat' using 56:57 with lines linecolor 4 notitle, \
'kadai3-2_sys3.dat' using 58:59 with lines linecolor 4 notitle, \
'kadai3-2_sys3.dat' using 60:61 with lines linecolor 4 notitle, \
'kadai3-2_sys3.dat' using 62:63 with lines linecolor 4 notitle, \
'kadai3-2_sys3.dat' using 64:65 with lines linecolor 4 notitle, \
'kadai3-2_sys3.dat' using 66:67 with lines linecolor 5 title 'r=1.50', \
'kadai3-2_sys3.dat' using 68:69 with lines linecolor 5 notitle, \
'kadai3-2_sys3.dat' using 70:71 with lines linecolor 5 notitle, \
'kadai3-2_sys3.dat' using 72:73 with lines linecolor 5 notitle, \
'kadai3-2_sys3.dat' using 74:75 with lines linecolor 5 notitle, \
'kadai3-2_sys3.dat' using 76:77 with lines linecolor 5 notitle, \
'kadai3-2_sys3.dat' using 78:79 with lines linecolor 5 notitle, \
'kadai3-2_sys3.dat' using 80:81 with lines linecolor 5 notitle
