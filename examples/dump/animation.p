set term gif animate delay 5 size 400, 400
set output "point.gif"
FN = "data_position.txt"
delta = 0.5
v(x,y) = sqrt(x**2 + y**2 + delta)

#set dgrid3d
#set hidden3d

set contour base
#unset surface
set view map
#set isosamples 4
set cntrparam levels incr -4,1,5
set xrange [-5:5]
set yrange [-5:5]

do for [n=1:400] {
    splot v(x,y) w l ls 1 nosurface, \
    FN u 1:2:(0) every ::::n w lp lw 1 nocontour t sprintf("n=%i", n)
}
