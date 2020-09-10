reset
set terminal postscript landscape enhanced color "Times_Roman" 32
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95

set output "fig01.eps"

funit=1.000000
set key left top
set xlabel "Phonon frequency (THz)"
set ylabel "Phonon DOS (1/THz)" 1
plot 'V1.000/vdos.out' using ($1*funit*1.e-12):($2/funit*1.e12) title "V/V0=1.000" w l lt -1 lw 1, \
 'V1.015/vdos.out' using ($1*funit*1.e-12):($2/funit*1.e12) title "V/V0=1.015" w l lt 1 lw 3

set output "fig02.eps"

set key left top
set xlabel "Band energy - E_F (eV)"
set ylabel "Electronic DOS (1/eV)" 1
plot [x=-10:10] 'V1.000/t.dat' using ($1-7.81774960):($2+$3) notitle w l lt -1 lw 1, \
 'vline.dat' using 1:2 title "E_F" w l lt 4 lw 3
