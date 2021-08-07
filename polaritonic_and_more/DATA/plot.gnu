#!/usr/local/bin/gnuplot
set terminal postscript eps enhanced color size 3in,3in 'Helvetica'
set output 'MgH_ccPVDZ_Ez_5mH.eps'

set size square
set encoding iso_8859_1
set key bottom right
#set title 'Potential Energy Surfaces of MgH+ Ez = 5 mH'
set ylabel 'Energy (Hartrees)'
set xlabel 'Bondlength ({\305})'
plot 'MgH_ccPVDZ_Ez_5mH.txt' u 1:2 w l lw 4 lc rgb "red" dt 1 title 'X,1', \
'MgH_ccPVDZ_Ez_5mH.txt' u 1:3 w l lw 4 lc rgb "blue" dt 1 title 'A,0', \
'MgH_ccPVDZ_Ez_5mH.txt' u 1:6 w p pt 7 lc rgb "red" title 'CQED-CIS LP', \
'MgH_ccPVDZ_Ez_5mH.txt' u 1:7 w p pt 7 lc rgb "blue" title 'CQED-CIS UP', \
'MgH_ccPVDZ_Ez_5mH.txt' u 1:4 w l lw 4 lc rgb 'red'  dt 2 title 'Model LP', \
'MgH_ccPVDZ_Ez_5mH.txt' u 1:5 w l lw 4 lc rgb "blue" dt 2 title 'Model UP'

set output 'MgH_ccPVDZ_Ez_1mH.eps'

plot 'MgH_ccPVDZ_Ez_1mH.txt' u 1:2 w l lw 4 lc rgb "red" dt 1 title 'X,1', \
'MgH_ccPVDZ_Ez_1mH.txt' u 1:3 w l lw 4 lc rgb "blue" dt 1 title 'A,0', \
'MgH_ccPVDZ_Ez_1mH.txt' u 1:6 w p pt 7 lc rgb "red" title 'CQED-CIS LP', \
'MgH_ccPVDZ_Ez_1mH.txt' u 1:7 w p pt 7 lc rgb "blue" title 'CQED-CIS UP', \
'MgH_ccPVDZ_Ez_1mH.txt' u 1:4 w l lw 4 lc rgb 'red'  dt 2 title 'Model LP', \
'MgH_ccPVDZ_Ez_1mH.txt' u 1:5 w l lw 4 lc rgb "blue" dt 2 title 'Model UP'
