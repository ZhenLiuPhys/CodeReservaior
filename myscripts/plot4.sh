#!/bin/sh
gnuplot << EOF
set terminal png enhanced fon Helvetica 12 size 1200,900
set output "param_ma$1.png"
set size 1,1
set origin 0,0
set multiplot

set size 0.33, 0.33
set origin 0, 0.66
set xlabel 'm_A'
set ylabel '{/Symbol l}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23 : NaN): 1 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): 1 ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0.33, 0.66
set xlabel 'm_A'
set ylabel '{/Symbol k}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23 : NaN): 2 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): 2 ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0.66, 0.66
set xlabel 'm_A'
set ylabel 'tan{/Symbol b}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23 : NaN): 3 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): 3 ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0, 0.33
set xlabel 'm_A'
set ylabel '{/Symbol m}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23: NaN): 4 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): 4 ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0.33, 0.33
set xlabel 'm_A'
set ylabel 'A_{/Symbol l}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23: NaN): 5 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): 5 ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0.66, 0.33
set xlabel 'm_A'
set ylabel 'A_{/Symbol k}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23: NaN): 6 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): 6 ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0, 0
set xlabel 'm_A'
set ylabel 'M_{SQ3}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23: NaN): (sqrt(\$7)) ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): (sqrt(\$7)) ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0.33, 0
set xlabel 'm_A'
set ylabel 'M_{SU3}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23: NaN): (sqrt(\$8)) ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): (sqrt(\$8)) ti 'h_1=H_{SM}'

set size 0.33, 0.33
set origin 0.66, 0
set xlabel 'm_A'
set ylabel 'A_{t}'
unset key
plot "nmssmout$1.dat" u (\$27 < 127 && \$27 > 123 ? \$23: NaN): 12 ti 'h_2=H_{SM}',"nmssmout$1.dat" u (\$26 < 127 && \$26 > 123 ? \$23 : NaN): 12 ti 'h_1=H_{SM}'
EOF

