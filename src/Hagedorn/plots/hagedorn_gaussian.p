fbuild(ix,delta) = sprintf("%d%.3f",ix, delta)
str_index(ix) = sprintf("%d",ix)
fout = '~/NAQMD/src/Hagedorn/plots/hagedorn_gaussian.tex'
set output fout
set terminal cairolatex pdf size 6,7


load "../../Auxiliary/settings.p"

TITLE = 'Hagedorn (transmitted)wavepacket'
ROWS = 4
COLS = 2
set multiplot layout ROWS, COLS title TITLE font ",12"

#-------------------- PLOT 1 ---------------------------------------
set format x ''
array DELTA[COLS] = [0.1, 0.5] 
array IX[ROWS] = [0, 1, 3, 10]
do for [i=1:ROWS]{
  do for [j=1:COLS]{
    plot \
      "~/NAQMD/src/Hagedorn/data/hag_wave".fbuild(IX[i], DELTA[j])."_transmitted.txt" \
      u 1:($2) with lines t '$\hat{\varphi}^+_{'.str_index(IX[i]).'}$', \
      "" u 1:(-$4) with lines t '$\hat{\varphi}^-_{'.str_index(IX[i]).'}$'
  }
  if (i == 3) {set format x '%.1t'} 
}
unset multiplot
set output 
unset terminal
