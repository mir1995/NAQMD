fbuild(ix,delta) = sprintf("%d%.3f",ix, delta)
fbuildnew(delta, qc, eps, ix) = sprintf("delta%.3fqc%.3feps%.3fp4Hagedorn%d", delta, qc, eps, ix)
str_index(ix) = sprintf("%d",ix)
fout = '~/NAQMD/src/Hagedorn/plots/hagedorn_transmitted_exact.tex'
set output fout
set terminal cairolatex pdf size 6,7

load "../../Auxiliary/settings.p"

TITLE = 'Hagedorn (transmitted)wavepacket vs exact solution'
ROWS = 4
COLS = 2
set multiplot layout ROWS, COLS title TITLE font ",12"

#-------------------- PLOT 1 ---------------------------------------
set format x ''
array DELTA[COLS] = [0.1, 0.5] 
array IX[ROWS] = [0, 1, 3, 10]
set xrange [2:6]
do for [i=1:ROWS]{
  do for [j=1:COLS]{
    plot \
      "~/QMD/src/wip_flatevls/data/wip_simple1".fbuildnew(DELTA[j], 0.5, 0.1, IX[i]).".txt" \
      u 2:($5) with lines t '$\hat{\varphi}^-_{'.str_index(IX[i]).'}$', \
      "~/NAQMD/src/Hagedorn/data/hag_wave".fbuild(IX[i], DELTA[j])."_transmitted.txt" \
      u 1:(-$4) with lines t '$\hat{\varphi}^+_{'.str_index(IX[i]).'}$'
  }
  if (i == 3) {set format x '%.1t'} 
}
unset multiplot
set output 
unset terminal
