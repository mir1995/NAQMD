fout = '~/NAQMD/src/SurfaceHopping/plots/mass.tex'
set output fout
set terminal cairolatex pdf size 6,5


load "../../Auxiliary/settings.p"

TITLE = 'Mass Absolute Error'


set size 1,1
set autoscale fix
#set colorbox vertical user origin .02, .1 size .04,.5
set colorbox horizontal user origin 0.2, 0.4 size 0.6, 0.03

set view map
set dgrid3d

set multiplot layout 2,2 title TITLE font ",12" \
#-------------------- PLOT 1 ---------------------------------------
#set view 30,220,2,3


array COLS[2] = [3,5]
set xlabel "$\epsilon$"
set ylabel "$\delta$"

splot \
  "~/NAQMD/src/SurfaceHopping/plots/mass.txt" using 1:2:(abs($3)) with pm3d t 'goddard'
unset colorbox

splot \
  "~/NAQMD/src/SurfaceHopping/plots/mass.txt" using 1:2:(abs($5)) with pm3d t 'lasser'
unset multiplot
set output 
unset terminal

