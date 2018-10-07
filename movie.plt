#
set terminal pngcairo font "arial,24" size 800, 600 enhanced rounded truecolor

set lmargin at screen 0.15
set rmargin at screen 0.85
set bmargin at screen 0.20
set tmargin at screen 0.85

# color definitions
set style line 11 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 5 # --- red
set style line 12 lc rgb '#FF4500' pt 6 ps 1 lt 1 lw 5 # --- orangered
set style line 13 lc rgb '#B22222' pt 6 ps 1 lt 1 lw 5 # --- firebrick
set style line 14 lc rgb '#DC143C' pt 6 ps 1 lt 1 lw 5 # --- crimson

set style line 21 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 5 # --- green
set style line 22 lc rgb '#006400' pt 6 ps 1 lt 1 lw 5 # --- darkgreen
set style line 23 lc rgb '#228B22' pt 6 ps 1 lt 1 lw 5 # --- forestgreen
set style line 24 lc rgb '#808000' pt 6 ps 1 lt 1 lw 5 # --- olive

set style line 31 lc rgb '#8A2BE2' pt 6 ps 1 lt 1 lw 5 # --- blueviolet
set style line 32 lc rgb '#00008B' pt 6 ps 1 lt 1 lw 5 # --- darkblue

set style line 41 lc rgb '#2F4F4F' pt 6 ps 1 lt 1 lw 5 # --- darkslategray

set border 31 lw 2

nsteps = 50
dt = 1

model = 3
cases = 0

do for [i=0:cases]  {
  do for [j=0:99] {
    set output 'res'.i.'/wfc.state0.frame'.j.'.png'
    set xrange[-20:20]
    set yrange[0:0.01]
    set ylabel "Amplitude (a.u)"  offset 0.0, 0.5
    set xlabel "Position (a.u), t = ".(j*dt*nsteps)." (a.u) "  offset 0.0, 0.2
    set nokey
    plot  'res'.i.'/wfc.state0.frame'.j.'' u ($1):($2*0.005) w l lt 8 lw 8,\
          "_pes.txt" u ($1):($2) w l ls 11 lw 5 t "",\
  }
}
  

#do for [i=0:1]  {
#  do for [j=0:cases]  {
#    do for [k=0:74] {
#      set output '_1D_dist_qc/_dist'.model.''.i.''.j.'_t'.(k*dt*nsteps).'.png'
#      set xrange[-5:5]
#      set yrange[0:0.02]
#      set ylabel "Potential Energy (a.u)"  offset 0.0, 0.5
#      set xlabel "Position (a.u), t = ".(k*dt*nsteps)." (a.u) "  offset 0.0, 0.2
#      set nokey
#      plot '_1D_dist_qc/_dist_'.model.''.i.''.j.'_'.(k*dt*nsteps).'.txt' u 1:($2*0.1) w l lt 8 lw 8 t "",\
# #           "_pes.txt" u ($1):($2) w l ls 11 lw 5 t "",\
#
#    }
#  }
#}
#




