#!/usr/bin/gnuplot

#set terminal png transparent nocrop enhanced size 1750,1250 font "arial,20" 
#set output 'BSabs20deg.png'
reset
set isosample 150, 150
set cntrparam order 8
set contour surface
set cntrparam level discrete   410000, 440000
set cntrlabel start 25 interval -10 font ",7"
set view 60, 30, 1, 1.1
unset surface

set table 'cont.dat'
filter(x)=(x>-400) ? (x) : (1/0)
splot  "BFEdgeATEST" u 1:2:(filter($3))
unset table


reset
#set xrange [-1.3:1.3]
#set yrange [-1.3:1.3]
#set xrange [-1:1]
#set yrange [-1:1]
set xrange [-0.1:0.1]
set yrange [-0.1:0.1]
unset key

set palette rgbformulae  34,35,36
set style line 1 linecolor rgb "#332288"
set style line 2 linecolor rgb "#117733"

set style increment user
#set pm3d map
#set pm3d explicit
#set style textbox opaque noborder  
set contour
set cntrlabel onecolor start 25 interval 20 font ",7"
set cbrange [6000:9000] 
l '<./cont.sh cont.dat 0 10 0'
p  'BFEdgeATEST' u 1:2:3 with ima, '<./cont.sh cont.dat 1 10 0' w l lt 3 lw 2


