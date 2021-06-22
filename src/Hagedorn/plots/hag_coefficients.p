fbuild(delta) = sprintf("%.3f",delta)
str_index(ix) = sprintf("%d",ix)
fout = '~/NAQMD/src/Hagedorn/plots/hag_coefficients.tex'
set output fout
set terminal cairolatex pdf size 6.5,3.5


load "../../Auxiliary/settings.p"

TITLE = 'Coefficients for transmitted wavepacket'
ROWS = 1
COLS = 2
set multiplot layout ROWS, COLS title TITLE font ",12" \
      margins 0.1,0.9,0.1,0.9 \
      spacing 0.03,0.08
#-------------------- PLOT 1 ---------------------------------------
set xlabel '$l$'
array DELTA[COLS] = [0.1, 0.5] 
do for [j=1:COLS]{
  if (j==2){set format y ''} else {set ylabel '$|c_l|$'}
  plot \
    "~/NAQMD/src/Hagedorn/data/hag_coefficients0".fbuild(DELTA[j])."_noshift.txt" \
    every ::1::6 u 1:(abs($2 + {0,1}*$3)) w lp ls 1 t '$p_1^\prime = p$', \
    "~/NAQMD/src/Hagedorn/data/hag_coefficients0".fbuild(DELTA[j])."_shift.txt" \
    every ::1::6 u 1:(abs($2 + {0,1}*$3)) w lp ls 2 t '$p_1^\prime = \langle \xi \hat{\psi}, \hat{\psi}\rangle$', \
    "~/NAQMD/src/Hagedorn/data/hag_coefficients0".fbuild(DELTA[j])."_conservation.txt" \
    every ::1::6 u 1:(abs($2 + {0,1}*$3)) w lp ls 3 t '$p_1^\prime = p_1 + \sqrt{p_1^2 + 4\delta}$'
}
unset multiplot
set output 
unset terminal
