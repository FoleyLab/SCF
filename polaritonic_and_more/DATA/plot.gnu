#!/usr/local/bin/gnuplot
set terminal postscript eps enhanced color size 3in,3in 'Helvetica'
set output 'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_10g.eps'

set size square
set encoding iso_8859_1
set key bottom right
#set title 'Potential Energy Surfaces of MgH+ Lamz = 5 mH'
set ylabel 'Energy (Hartrees)'
set xlabel 'Bondlength ({\305})'
plot 'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_10g.txt' u 1:2 w l lw 4 lc rgb "red" dt 1 title 'X,1', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_10g.txt' u 1:3 w l lw 4 lc rgb "blue" dt 1 title 'A,0', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_10g.txt' u 1:8 w p pt 7 lc rgb "red" title 'CQED-CIS LP', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_10g.txt' u 1:9 w p pt 7 lc rgb "blue" title 'CQED-CIS UP', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_10g.txt' u 1:6 w l lw 4 lc rgb 'red'  dt 2 title 'Model LP', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_10g.txt' u 1:7 w l lw 4 lc rgb "blue" dt 2 title 'Model UP'


set output 'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.eps'

set size square
set encoding iso_8859_1
set key bottom right
#set title 'Potential Energy Surfaces of MgH+ Lamz = 5 mH'
set ylabel 'Energy (Hartrees)'
set xlabel 'Bondlength ({\305})'
plot 'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:2 w l lw 4 lc rgb "red" dt 1 title 'X,1', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:3 w l lw 4 lc rgb "blue" dt 1 title 'A,0', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:8 w p pt 7 lc rgb "red" title 'CQED-CIS LP', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:9 w p pt 7 lc rgb "blue" title 'CQED-CIS UP', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:6 w l lw 4 lc rgb 'red'  dt 2 title 'Model LP', \
'MgH_ccpVDZ_Lamz_12.5mH_om_4.75_0g.txt' u 1:7 w l lw 4 lc rgb "blue" dt 2 title 'Model UP'

set output 'MgH_gs_lamz_7.5mH_om_4.75.eps'
set size square
set encoding iso_8859_1
set key bottom right
#set title 'Potential Energy Surfaces of MgH+ Lamz = 5 mH'
set xrange [1.25:2.75]
set ylabel 'Energy (Hartrees)'
set xlabel 'Bondlength ({\305})'
plot 'MgH_gs_lamz_7.5mH_om_4.75.txt' u 1:2 w l lw 4 lc rgb "red" dt 1 title 'RHF', \
'MgH_gs_lamz_7.5mH_om_4.75.txt' u 1:3 w l lw 4 lc rgb "blue" dt 1 title 'CQED-RHF', \
'MgH_gs_lamz_7.5mH_om_4.75.txt' u 1:4 w l lw 4 lc rgb "dark-green" dt 1 title 'CQED-CIS'


reset
set terminal postscript eps enhanced color size 3in,3in 'Helvetica'
set output 'Formaldehyde_total_E.eps'
set size square
set encoding iso_8859_1
set key bottom right
set ylabel 'Energy (Hartrees)'
set xlabel 'Electric Field Strength (a.u.)'
plot 'Formaldehyde_ccpVDZ_variable_L.txt' u 1:2 w l lw 4 lc rgb "red" dt 1 title "E_y", \
'Formaldehyde_ccpVDZ_variable_L.txt' u 1:10 w l lw 4 lc rgb "blue" dt 1 title "E_z", \
'Formaldehyde_ccpVDZ_variable_L.txt' u 1:18 w l lw 4 lc rgb "dark-green" dt 1 title "E_y, E_z"

set terminal postscript eps enhanced color size 3in,3in 'Helvetica'
set output 'Formaldehyde_total_RHF_E.eps'
set size square
set encoding iso_8859_1
set key bottom right
set ylabel 'Energy (Hartrees)'
set xlabel 'Electric Field Strength (a.u.)'
plot 'Formaldehyde_ccpVDZ_variable_L.txt' u 1:($3+$4) w l lw 4 lc rgb "red" dt 1 title "E_y", \
'Formaldehyde_ccpVDZ_variable_L.txt' u 1:($11+$12) w l lw 4 lc rgb "blue" dt 1 title "E_z", \
'Formaldehyde_ccpVDZ_variable_L.txt' u 1:($19+$20) w l lw 4 lc rgb "dark-green" dt 1 title "E_y, E_z"

set terminal postscript eps enhanced color size 3in,3in 'Helvetica'
set output 'Formaldehyde_total_PF_E.eps'
set size square
set encoding iso_8859_1
set key bottom right
set ylabel 'Energy (Hartrees)'
set xlabel 'Electric Field Strength (a.u.)'
plot 'Formaldehyde_ccpVDZ_variable_L.txt' u 1:($5+$6+$7) w l lw 4 lc rgb "red" dt 1 title "E_y", \
'Formaldehyde_ccpVDZ_variable_L.txt' u 1:($13+$14+$15) w l lw 4 lc rgb "blue" dt 1 title "E_z", \
'Formaldehyde_ccpVDZ_variable_L.txt' u 1:($21+$22+$23) w l lw 4 lc rgb "dark-green" dt 1 title "E_y, E_z"



#plt.plot(l_mag, cqed_rhf_energy_array[1,:]+cqed_rhf_energy_array[2,:], label="E_y" )
#plt.plot(l_mag, cqed_rhf_energy_array[9,:]+cqed_rhf_energy_array[10,:], label="E_z" )
#plt.plot(l_mag, cqed_rhf_energy_array[17,:]+cqed_rhf_energy_array[18,:], label="E_y, E_z" )
#plt.legend()
#plt.show()

#plt.plot(l_mag, cqed_rhf_energy_array[3,:]+cqed_rhf_energy_array[4,:]+cqed_rhf_energy_array[5,:], label="E_y" )
#plt.plot(l_mag, cqed_rhf_energy_array[11,:]+cqed_rhf_energy_array[12,:]+cqed_rhf_energy_array[13,:], label="E_z" )
#plt.plot(l_mag, cqed_rhf_energy_array[19,:]+cqed_rhf_energy_array[20,:]+cqed_rhf_energy_array[21,:], label="E_y, E_z" )
#plt.legend()
#plt.show()

#plt.plot(l_mag, cqed_rhf_energy_array[0,:], label="E_y" )
#plt.plot(l_mag, cqed_rhf_energy_array[8,:], label="E_z" )
#plt.plot(l_mag, cqed_rhf_energy_array[16,:], label="E_y, E_z" )
