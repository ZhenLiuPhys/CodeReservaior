#!/bin/sh
gnuplot << EOF
set terminal png enhanced fon Helvetica 12 size 800,600
set output "couplings$1.png"
set size 1,1
set origin 0,0
set multiplot

set size 0.5, 0.5
set origin 0, 0.5
set xlabel 'H_{SM}WW'
set ylabel 'H_{SM}uu'
set key top center
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$88 : NaN): 86 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$83 : NaN): 81 ti 'h_1=H_{SM}'

set size 0.5, 0.5
set origin 0.5, 0.5
set xlabel 'H_{SM}WW'
set ylabel 'H_{SM}dd'
set key top center
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$88 : NaN): 87 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$83 : NaN): 82 ti 'h_1=H_{SM}'

set size 0.5, 0.5
set origin 0, 0
set xlabel 'H_{SM}WW'
set ylabel 'H_{SM}gg'
set key top center
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$88 : NaN): 89 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$84 : NaN): 82 ti 'h_1=H_{SM}'

set size 0.5, 0.5
set origin 0.5, 0
set xlabel 'H_{SM}WW'
set ylabel 'H_{SM}{/Symbol gg}'
set key top center
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$88 : NaN): 90 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$85 : NaN): 82 ti 'h_1=H_{SM}'
EOF

