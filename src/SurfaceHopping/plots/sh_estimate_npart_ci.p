fbuild(p) = sprintf("%d",p)
fout = '~/NAQMD/src/SurfaceHopping/plots/sh_estimate_npart_ci.tex'

set output fout
set terminal cairolatex pdf size 4,3

load "../../Auxiliary/settings.p"

TITLE = 'npart with 95'

#-------------------- PLOT 1 ---------------------------------------
set xlabel '$\epsilon$'
set ylabel '$\sigma\left(a \circ \Phi^{t_f}\right)$'

set logscale xy
set key horizontal left top

plot for [p=1:6] \
  "~/NAQMD/src/SurfaceHopping/data/variancedelta0.5000alpha0.5000npart10000.txt" \
  u 1:($2==p?$3:1/0) t '$p=$'.fbuild(p)

set output 
unset terminal
