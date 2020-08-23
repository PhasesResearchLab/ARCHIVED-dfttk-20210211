reset
set terminal postscript landscape enhanced color "Times_Roman" 20
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95

set output "vdis.eps"

qunit=1.0
eunit=1.000000
funit=1.000000
p0 = 0.000000
pp0 = 0.183232
p1 = 0.183232
pp1 = 0.129565
p2 = 0.312797
pp2 = 0.129565
p3 = 0.442361
pp3 = 0.107512
p4 = 0.549874
pp4 = 0.129565
p5 = 0.679438
pp5 = 0.129565
p6 = 0.809003
pp6 = 0.212445
p7 = 1.021448

set xtics ( '{/Symbol G}' qunit*0.000000, 'M' qunit*0.183232, 'X' qunit*0.312797, '{/Symbol G}' qunit*0.442361, 'Z' qunit*0.549874, 'A' qunit*0.679438, 'R' qunit*0.809003, '{/Symbol G}' qunit*1.021448)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p7*1.0001] [funit*-2.000000:funit*25.000000] \
'vline.dat' using (qunit*$1):(funit*$2) notitle w l lt 4, \
 'vdis.out' index 0 using (qunit*p0+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 4 using (qunit*p4+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 5 using (qunit*p5+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 6 using (qunit*p6+qunit*$1):(funit*$19) notitle w l lt -1

#
#
