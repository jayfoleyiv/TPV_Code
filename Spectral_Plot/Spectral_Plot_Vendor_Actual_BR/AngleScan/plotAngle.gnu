#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'TE_vs_Angle.eps'

set xlabel 'Wavelength (nm)'
set ylabel 'Thermal Emission'
plot '0Deg_BR.txt' u 1:6 w l lw 2 title 'Normal', \
'5Deg_BR.txt' u 1:6 w l lw 2 title '5 Deg', \
'15Deg_BR.txt' u 1:6 w l lw 2 title '15 Deg', \
'25Deg_BR.txt' u 1:6 w l lw 2 title '25 Deg', \
'35Deg_BR.txt' u 1:6 w l lw 2 title '35 Deg', \
'45Deg_BR.txt' u 1:6 w l lw 2 title '45 Deg', \
'0Deg_BR.txt' u 1:7 w l lt -1 title 'Blackbody, 1700 K'


