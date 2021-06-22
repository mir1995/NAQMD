fbuild(p) = sprintf("%d",p)
fout = '~/NAQMD/src/SurfaceHopping/plots/sh_trajectories.tex'

set output fout
set terminal cairolatex pdf size 6,9

load "../../Auxiliary/settings.p"

TITLE = 'Absolute error - Mass, Energy, $\delta \in \{0.05, 0.5\}$'

#-------------------- PLOT 1 ---------------------------------------

set logscale y
set key horizontal left top

set multiplot layout 5,2 title TITLE font ",12" 

set lmargin 1.5
set bmargin 1.5
set rmargin 0
set tmargin 1

set xtics format " "
set ylabel '$\|\phi^+\|_{L^2}$'
set format y "10^{%L}"

set title 'Mass Up $\delta = 0.05$'

set key vertical

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($10 - $21)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($10 - $21)) t 'lasser'

set ytics format " "
unset ylabel
unset key

set title 'Mass Down $\delta = 0.5$'

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($10 - $21)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($10 - $21)) t 'lasser'

set ytics
set format y "10^{%L}"

#############################################################
set title 'Mean Position Up $\delta = 0.05$'

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($2 - $13)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($2 - $13)) t 'lasser'

set ytics format " "
set title 'Mean Position Up $\delta = 0.5$'

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($2 - $13)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($2 - $13)) t 'lasser'

#############################################################
set title 'Mean Position Down $\delta = 0.05$'
set ytics
set format y "10^{%L}"

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($3 - $14)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($3 - $14)) t 'lasser'

set ytics format " "
set title 'Mean Position Down $\delta = 0.5$'

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($3 - $14)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($3 - $14)) t 'lasser'
set xtics
set xlabel 'itr'
set ytics
set format y "10^{%L}"
#############################################################
set title 'Energy Down $\delta = 0.05$'

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($9 - $20)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.0500alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.0500eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($9 - $20)) t 'lasser'

set ytics format " "

set title 'Energy Down $\delta = 0.5$'

plot \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000goddard.txt" \
       u 1:(abs($9 - $20)) t 'goddard', \
       "< paste ~/QMD/src/wip_simple1/data/observablesdelta0.5000alpha0.5000eps0.1000p4.txt ~/NAQMD/src/SurfaceHopping/data/observablesdelta0.5000eps0.1000alpha0.5000p4fnpart100000lasser.txt" \
       u 1:(abs($9 - $20)) t 'lasser'


unset multiplot
set output 
unset terminal
