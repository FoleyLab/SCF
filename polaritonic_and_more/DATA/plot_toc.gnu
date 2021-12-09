#!/usr/local/bin/gnuplot
set terminal postscript eps enhanced color 
set output 'test.eps'

set size square
set encoding iso_8859_1
unset key
unset border
unset xtics
unset ytics
#set title 'Potential Energy Surfaces of MgH+ Lamz = 5 mH'
plot 'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:8 w l smooth csplines lw 7 lc rgb "goldenrod" notitle, \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:9 w l smooth csplines lw 7 lc rgb "dark-orange" notitle

