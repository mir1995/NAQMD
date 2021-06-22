fout = '~/NAQMD/src/SurfaceHopping/plots/energy.tex'
set output fout
set terminal cairolatex pdf size 6.5,3.5


load "../../Auxiliary/settings.p"

TITLE = 'Energy Absolute Error'

set multiplot layout 2,2 title TITLE font ",12" \
      margins 0.1,0.9,0.1,0.9 \
      spacing 0.03,0.08
#-------------------- PLOT 1 ---------------------------------------
set view map
set dgrid3d
do for [j=1:4]{
  if (j==2){set format y ''} else {set ylabel '$|c_l|$'}
  splot \
    "~/NAQMD/src/SurfaceHopping/plots/energy.txt" using 1:2:3 with pm3d t ''
}
unset multiplot
set output 
unset terminal

