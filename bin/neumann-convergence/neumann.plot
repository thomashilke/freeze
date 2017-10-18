
set logscale xy

set xlabel 'h'
set ylabel 'Error L2'

plot 'neumann.dat' u (1./$1):2 w lp title 'Error L2', \
     '' u (1./$1):(1.0 * (1./$1)) w l title 'Slope 1', \
     '' u (1./$1):(0.1 * (1./$1)**(1./2.)) w l title 'Slope 1/2', \
     '' u (1./$1):(((100.)/(100.**(4./3.))/$1)**(3./4.)) w l title 'Slope 3/4'