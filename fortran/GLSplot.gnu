set term post por 14 size 22 cm, 15 cm solid color; set out "GLS.ps";
set macro
set border 31 lw 0

# plot frequency spectrum
load "GLS.log"
set title "Keplerian Periodogram for ".file offset -1
@Kep set title "Periodograms and windowfunction for ".file offset -1
set xlabel "Frequency [1/d]"; set ylabel "Power p"
set mxtics; set mytics
set xzeroaxis lt -1 lw 0

plkep=""; @kep plkep='"GLSKep.plt"  t "" w l lt 1 @Kep, '
plot @plkep "GLS.plt" t "" w l lt 3, "" us 1:(($3-1.)/5.) t "" w l lt 2

set out
set term x11  0 ti "GLS" font "arial,15";
repl
#reset

!epstopdf GLS.ps
#!echo "kpdf GLS.pdf &";

