# set terminal png transparent nocrop enhanced size 450,320 font "arial,8" 
# set output 'contours.23.png'
set border 15 front lt black linewidth 1.000 dashtype solid
set key at screen 1, 0.9, 0 right top vertical Right noreverse enhanced autotitle nobox
unset key
set style textbox opaque margins  0.5,  0.5 noborder
set logscale z 10
set view map scale 1
set isosamples 50, 50
unset surface 
set contour base
set cntrlabel onecolor format '%8.3g' font ',7' start 25 interval -1
set hidden3d back offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set cntrparam order 8
set cntrparam bspline
set cntrparam levels discrete 0.1,1 ,10 ,100 
set style data lines


set xr [-2.5:2.5]
set yr [-2.5:2.5]
unset xlabel
unset ylabel


set multiplot

set isosamples 500,100
filter(x)=(x>-400) ? (x) : (1/0)
splot  "contrast_edgeBd_fixed" u 1:2:(filter($3)) with lines lc rgb "#007700"

set isosamples 50,50
set cntrlabel start 25 interval -1 font ",7"
l '<./cont.sh cont.dat 0 13 0'
p  'contrast_edgeBd_fixed' u 1:2:3 with ima, '<./cont.sh cont.dat 1 13 0' w l lt 3 lw 2

unset multiplot
