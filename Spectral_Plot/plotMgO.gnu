#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25

set xlabel 'Wavelength (nm)'
set ylabel 'Thermal Emission'

set output 'MgO_Cr2O3.eps'

plot 'MgO_Cr2O3_Paretto.txt' u 1:6 w l lw 2 title 'SE: 60, SD: 1.72e5', \
'MgO_Cr2O3_Paretto.txt' u 1:7 w l lw 2 title 'Blackbody', \
