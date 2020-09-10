reset
set terminal postscript landscape enhanced color "Times_Roman" 32
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95


set output "fig01.eps"

set key right bottom
set ylabel "LEC (10**-6 1/K)"
set xlabel "T (K)"
set key right bottom
set label 1 "Ni" at graph 0.06,0.86 left font "Times_Roman,42"
plot [x=0:1750] [0:30] \
     'vdos_e' using 1:6 title "Phonon" w l lt 1 lw 2, \
      'exp02-1.dat' index 0 using 1:3 title "Touloukian'75" w p pt 6, \
      'exp02-1.dat' index 1 using 1:2 title "Touloukian'89" w p pt 12

set output "fig02.eps"

set key right bottom
set ylabel "Heat capacities (J/mol)"
set xlabel "T (K)"
set label 1 "Ni" at graph 0.06,0.86 left font "Times_Roman,42"
plot [x=0:1750] [0:45] \
     'vdos_e' using 1:7 title "Phonon" w l lt 1 lw 2, \
      'exp02-2.dat' index 1 using 1:2 title "Desai" w p pt 6, \
      'exp06.dat' index 0 using 1:3 title "Barin" w p pt 7

set output
