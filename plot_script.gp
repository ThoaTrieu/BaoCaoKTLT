set title 'Do thi ham so f(x)'
set xlabel 'x'
set ylabel 'f(x)'
set grid
plot 'data.dat' using 1:2 with lines title 'f(x)'
pause -1
