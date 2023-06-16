#!/bin/sh
gnuplot << EOF
set terminal png enhanced fon Helvetica 12 size 640, 480
set style data dots
set output "masses$1.png"
set xlabel 'm_A (GeV)'
set ylabel 'm (GeV)'
set key top left
plot "nmssmout$1.dat" u 23:26 ti 'm_{h_1}(Red)', "nmssmout$1.dat" u 23:27 ti 'm_{h_2}(Green)', "nmssmout$1.dat" u 23:28 ti 'm_{h_3}(Blue)', "nmssmout$1.dat" u 23:29 ti 'm_{a_1}(Purple)',"nmssmout$1.dat" u 23:30 ti 'm_{a_2}(Lightblue)', "nmssmout$1.dat" u 23:31 ti 'm_{h^+}(Brown)'
EOF

gnuplot << EOF
set terminal png enhanced fon Helvetica 12 size 640, 480
set output "masses_z$1.png"
set xlabel 'm_A (GeV)'
set ylabel 'm (GeV)'
set xrange[0:1500]
set key top left
plot "nmssmout$1.dat" u 23:26 ti 'm_{h_1}(Red)', "nmssmout$1.dat" u 23:27 ti 'm_{h_2}(Green)', "nmssmout$1.dat" u 23:28 ti 'm_{h_3}(Blue)', "nmssmout$1.dat" u 23:29 ti 'm_{a_1}(Purple)',"nmssmout$1.dat" u 23:30 ti 'm_{a_2}(Lightblue)', "nmssmout$1.dat" u 23:31 ti 'm_{h^+}(Brown)'
EOF
