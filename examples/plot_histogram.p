n=100 #number of intervals
max=4. #max value
min=-4. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set term png #output terminal and file
set output "histogram.png"
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "x"
set ylabel "Relative Frequency"
#count and plot
sum = 0
s(x) = ((sum = sum + 1), 0)
plot "position.txt" u ($1):(s($1))
plot "position.txt" u (hist($1,width)):(1.0/(width*sum)) smooth freq w boxes lc rgb"green" notitle
set output
