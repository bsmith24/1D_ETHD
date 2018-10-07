#
set terminal pngcairo font "arial,24" size 800, 600 enhanced rounded truecolor

set lmargin at screen 0.22
set rmargin at screen 0.95
set bmargin at screen 0.16
set tmargin at screen 0.95

# color definitions
set style line 11 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 5 # --- red
set style line 12 lc rgb '#FF0000' pt 1 ps 1 lt 1 lw 5 # --- bright red
set style line 13 lc rgb '#FF4500' pt 6 ps 1 lt 1 lw 5 # --- orangered
set style line 14 lc rgb '#B22222' pt 6 ps 1 lt 1 lw 5 # --- firebrick
set style line 15 lc rgb '#DC143C' pt 6 ps 1 lt 1 lw 5 # --- crimson
set style line 16 lc rgb '#FF7000' pt 6 ps 1 lt 1 lw 5 # --- dark orange

set style line 21 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 5 # --- green
set style line 22 lc rgb '#006400' pt 6 ps 1 lt 1 lw 5 # --- darkgreen
set style line 23 lc rgb '#228B22' pt 6 ps 1 lt 1 lw 5 # --- forestgreen
set style line 24 lc rgb '#808000' pt 6 ps 1 lt 1 lw 5 # --- olive
set style line 25 lc rgb '#00FF10' pt 6 ps 1 lt 1 lw 5 # --- lime
set style line 26 lc rgb '#20B2AA' pt 6 ps 1 lt 1 lw 5 # --- light sea green

set style line 31 lc rgb '#8A2BE2' pt 6 ps 1 lt 1 lw 5 # --- blueviolet
set style line 32 lc rgb '#000E8B' pt 6 ps 1 lt 1 lw 5 # --- royalblue
set style line 33 lc rgb '#00008B' pt 6 ps 1 lt 1 lw 5 # --- darkblue
set style line 34 lc rgb '#800080' pt 6 ps 1 lt 1 lw 5 # --- purple

set style line 41 lc rgb '#2F4F4F' pt 6 ps 1 lt 1 lw 5 # --- darkslategray
set style line 42 lc rgb '#C0C0C0' pt 6 ps 1 lt 1 lw 5 # --- silver
set style line 43 lc rgb '#D2691E' pt 6 ps 1 lt 1 lw 5 # --- chocolate



set linetype 1 dashtype (7,7,4,4,6,6) lc rgb '#FF0000' 
set linetype 2 dashtype (7,7,4,4,6,6) lc rgb '#5e9c36'
set linetype 3 dashtype (7,7,4,4,6,6) lc rgb '#000E8B'
set linetype 4 dashtype (7,7,4,4,6,6) lc rgb '#FF7000'
set linetype 5 dashtype (7,7,4,4,6,6) lc rgb '#800080'

set border 31 lw 2
set key font ",20"

model = 2
cases = 0

do for [i=0:cases]  {

unset xrange
unset yrange
set output "_fig2_".i."_model".model.".png"
set key top left
#set xrange[0:9]
#set yrange[-0.1:1.0]
set xlabel "t (fs)" offset 0.0, 0.5
set ylabel "<q>" offset 1.0, 0.0
plot   "_output_20".i.".txt" u ($1/41.3410):($6) w l lt 3 lw 6 t "",\
       "_output_21".i.".txt" u ($1/41.3410):($6) pt 7 ps 1 lc rgb '#000E8B' t "",\
       "_pops".i.".txt" u   ($1/41.3413745):($6) w l ls 32 lw 6 t  "",\

unset xrange
unset yrange
set output "_fig3_".i."_model".model.".png"
set key top left
#set xrange[0:9]
#set yrange[-0.1:1.0]
set xlabel "t (fs)" offset 0.0, 0.5
set ylabel "<p>" offset 1.0, 0.0
plot   "_output_20".i.".txt" u ($1/41.3410):($7) w l lt 3 lw 6 t "",\
       "_output_21".i.".txt" u ($1/41.3410):($7) pt 7 ps 1 lc rgb '#000E8B' t "",\
       "_pops".i.".txt" u   ($1/41.3413745):($7) w l ls 32 lw 6 t  "",\

unset xrange
unset yrange
set output "_fig4_".i."_model".model.".png"
set key top left
#set xrange[0:9]
#set yrange[-0.1:1.0]
set xlabel "t (fs)" offset 0.0, 0.5
set ylabel "<q^2>" offset 0.5, 0.0
plot   "_output_20".i.".txt" u ($1/41.3410):($8) w l lt 3 lw 6 t "",\
       "_output_21".i.".txt" u ($1/41.3410):($8) pt 7 ps 1 lc rgb '#000E8B' t "",\
       "_pops".i.".txt" u   ($1/41.3413745):($8) w l ls 32 lw 6 t  "",\

unset xrange
unset yrange
set output "_fig5_".i."_model".model.".png"
set key top left
#set xrange[0:9]
#set yrange[-0.1:1.0]
set xlabel "t (fs)" offset 0.0, 0.5
set ylabel "<p^2>" offset 0.5, 0.0
plot   "_output_20".i.".txt" u ($1/41.3410):($9) w l lt 3 lw 6 t "",\
       "_output_21".i.".txt" u ($1/41.3410):($9) pt 7 ps 1 lc rgb '#000E8B' t "",\
       "_pops".i.".txt" u   ($1/41.3413745):($9) w l ls 32 lw 6 t  "",\

}




