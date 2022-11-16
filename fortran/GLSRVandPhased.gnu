set nomultiplot
set size 1,1

set macro
set term x11  0 ti "GLS RV" font "arial,15"; set style line 1 pt 7 lt 3
# set term post por 14 size 22 cm, 15 cm solid color; set out "GLSRV.ps"; set style line 1 pt 7 lt 7;

# read last best-fitting parameters
load "GLS.log"
@igls  load "< sed 's@+/-.*$@;@'  iGLS[SK][ie][np].log"


# define Keplerian RV curve (parametric form)
M(E,e,T0,P) = (E-e*sin(E))/2/pi*P+T0
RV(E,RV0,K,e,w) = RV0+K*(e*cos(w*pi/180)+cos(w*pi/180+2*atan(((1+e)/(1-e))**0.5*tan(E/2))))

      gls="";  P=PSin; resfile="RVSinRes.plt";  fitcol="$2"; phasecol="$1";titex=""
@gls  RVfunc = "Amp*sin((x-Tph0)/P*2*pi)+Coff"
@kep  RVfunc = "M(E,e,T0,P), RV(E,RV0,K,e,w)"
@kep  gls="#"; P=PKep; resfile="RVKepRes.plt";  fitcol="$3"

@gls  @igls          resfile="iRVSinRes.plt"; fitcol="Amp*sin($1*2*pi)+Coff";
@kep  @igls          resfile="iRVKepRes.plt"; fitcol="RV(2*pi*$1,RV0,K,e,w)";  phasecol="(M(2*pi*$1,e,T0,P)-T0)/P"
@kep  set par;  set dummy E; titex=sprintf(", e=%.2f",e)

set mxtics; set mytics
err=""
@chi err=":3 w e"

# preparing multiplot
set multiplot
r=0.8; q=0.7; lmar0=9; tmar0=3

#### upper right plot #### RV phased ############
set size 1-r,q; set ori r,1-q; set tmar tmar0; set bmar 0; set lmar 0
set samp 200
unset xlabel; set form x '';
set xra [0:1]; set xtics 0,.5  nomirr; set mxtics 5
unset ylabel; set form y ''
set x2label "Phase [d]"; set x2ra[0:P]; set x2tics auto; set mx2tics
plot 'GLSFit.plt' us (@phasecol):(@fitcol)  w l t '',\
     resfile us 4:5@err t ""  ls 1
unset x2label

#### lower right plot #### residuals phased ###########
set size 1-r,1-q; set ori r,0; set tmar 0; set bmar 3; set lmar 0
set xlabel 'Phase'; set xtics mirr; set form x
set ylabel;
unset x2label; unset x2tics
set xzeroaxis lt -1 lw 0
plot resfile us 4:2@err t ""  ls 1


#### lower left plot #### residuals ###########
set size r,1-q; set ori 0,0; set tmar 0; set bmar 3; set lmar lmar0
set xlabel 'Time'; set xrange [*:*] ; set xtics auto; set form x; set form y
set ylabel 'O-C'
set xzeroaxis lt -1 lw 0
plot  resfile us 1:2@err t ""  ls 1

set xrange [GPVAL_X_MIN:GPVAL_X_MAX]
set yrange [*:*]
@kep set trange [(GPVAL_X_MIN-T0)/P*2*pi-pi: (GPVAL_X_MAX-T0)/P*2*pi+pi ]
cyc=(GPVAL_X_MAX-GPVAL_X_MIN)/P  # number of cycles
set samp 100*(cyc>1?cyc:1)


#### upper left plot #### RV time series ############
set size r,q; set ori 0,1-q; set tmar tmar0; set bmar 0; set lmar lmar0
unset xlabel; set form x "";
set ylabel 'RV'
unset xzeroaxis
plot  @RVfunc t sprintf("P=%f d",P).titex,\
      resfile us 1:5@err t ""  ls 1


unset multiplot
#set out
set form x; unset rmar; unset bmar; unset tmar; unset lmar; set xlabel 'Time';
#reset

unset size; unset ori; unset form;

#set term x11

#!echo "kpdf GLSRV.pdf &"

