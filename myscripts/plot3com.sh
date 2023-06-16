#!/bin/sh
gnuplot << EOF
set terminal png enhanced fon Helvetica 12 size 800,600
set output "param$1.png"
set size 1,1
set origin 0,0
set style data dots
set multiplot

set size 0.5, 0.5
set origin 0, 0.5
set xlabel '{/Symbol l}'
set ylabel '{/Symbol k}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$1 : NaN): 2 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$1 : NaN): 2 ti 'h_1=H_{SM}'

set size 0.5, 0.5
set origin 0.5, 0.5
set xlabel 'tan{/Symbol b}'
set ylabel '{/Symbol m}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$3 : NaN): 4 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$3 : NaN): 4 ti 'h_1=H_{SM}'

set size 0.5, 0.5
set origin 0, 0
set xlabel 'A_{/Symbol l}'
set ylabel 'A_{/Symbol k}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$5 : NaN): 6 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$5 : NaN): 6 ti 'h_1=H_{SM}'

set size 0.5, 0.5
set origin 0.5, 0
set xlabel 'M_{SQ_3}'
set ylabel 'M_{SU_3}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? (sqrt(\$7)) : NaN): (sqrt(\$8)) ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? (sqrt(\$7)) : NaN): (sqrt(\$8)) ti 'h_1=H_{SM}'
EOF

