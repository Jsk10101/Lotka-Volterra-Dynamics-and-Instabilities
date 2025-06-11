# This is our plotting file for the lotka_volterra.cpp file and the process that goes with running it. Any information need about describing this file is within the Overview.cpp file

set grid
set xlabel "Time"
set ylabel "Population" 
set title "Lotka Volterra Model"

plot "lotka_volterra.dat" u 1:2 w l title "Prey","lotka_volterra.dat" u 1:3 w l title "Predator"

replot
