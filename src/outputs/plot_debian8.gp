set terminal png size 1024, 768
set output 'newton_basins.png'
unset key
set style data lines
unset cbtics
unset ytics
unset xtics
set cbrange [-500:500] noreverse nowriteback # defines the range of point values and the mapping to the color
set palette rgbformulae 33, 13, 10

gnuplot --persist -e 'plot "output.txt" using 2:1:3 with image'
