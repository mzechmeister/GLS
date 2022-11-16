set nomultiplot
set macro
save term "plot.gnu.tmp2";
`awk '/x11/ {print "#"}' plot.gnu.tmp2` set term x11  0 ti "GLS" font "arial,15"; set style line 1 pt 7 lt 3 lc rgb "black"
`awk '/post/{print "#"}' plot.gnu.tmp2` set term post por 18 size 22 cm, 15 cm solid color; set out "GLSRV.ps"; set style line 1 pt 7 lt 7;

### load parameters
impr=""
load "GLS.log"
@impr `awk '$0=$1$2$3";"' GLSfit.log`;

### define Kepler  Keplerian RV curve
M(E,e,T0,P) = (E-e*sin(E))/2/pi*P+T0
RV(E,RV0,K,e,w) = RV0+K*(e*cos(w*pi/180)+cos(w*pi/180+2*atan(((1+e)/(1-e))**0.5*tan(E/2))))
RVcurve='M(E,e,T0,PKep)-JD0, RV(E,RV0,K,e,w)'

# define Kepler equation
Ex(x)=Ex_sol(x,x,x+e*sin(x)+e**2/2*sin(2*x))
# solve recursive
Ex_sol(x,F,Fn)=(abs((Fn-F)/F)<0.00001)? Fn:Ex_sol(x,Fn,Fn+(x-Fn+e*sin(Fn))/(1-e*cos(Fn))) ;

@impr RV(x,RV0,K,e,w,t0,P) = RV0 + K*(e*cos(w*pi/180) +cos(w*pi/180+2*atan(((1+e)/(1-e))**0.5*tan(Ex((x-t0)/P*2*pi)/2))))
@impr RVcurve='RV(x,RV0,K,e,w,T0,P)'
     gls="";  P=PSin; resfile="RVSinRes.plt"
@kep gls="#"; P=PKep; resfile="RVKepRes.plt"
@kep set par;  set dummy E;
@impr unset par
rescol=2;
@impr rescol='($5-Amp*sin(($1-Tph0)/PSin*2*pi)-Coff)'
@impr @kep rescol='($5-RV($1,RV0,K,e,w,T0,PKep))'
set mxtics; set mytics
set lmar 10

set multiplot
set size 1,0.3; set ori 0,0
set xlabel 'Time'
set ylabel 'O-C'
set form x
set tmar 0
set xzeroaxis lt -1 lw 0

### lower panel #### residuals ####
pl  resfile us 1:@rescol:3 w e t ""  ls 1

set xrange [GPVAL_X_MIN:GPVAL_X_MAX]
set yrange [*:*]
set trange [(GPVAL_X_MIN-t0)/P*2*pi-pi: (GPVAL_X_MAX-t0)/P*2*pi+pi ]
cyc=(GPVAL_X_MAX-GPVAL_X_MIN)/P
set samp 100*(cyc>1?cyc:1)

set size 1,0.7; set ori 0,0.3
set ylabel 'RV'
set form x ""
unset tmar
unset xlabel; unset xzeroaxis

### upper panel #### RV time series ####
@kep pl  @RVcurve lc rgb "red" tit sprintf("P=%f d, e=%.2f",PKep,e),\
@kep     resfile us 1:5:3 w e t ""  ls 1
@gls pl  Amp*sin((x-Tph0)/PSin*2*pi)+Coff    lc rgb "red" tit sprintf("P=%f d",PSin),\
@gls     resfile us 1:5:3 w e t ""  ls 1


unset multiplot
set out
#reset
`awk '/post/{print "#"}' plot.gnu.tmp2` ; reread   # read again for x11-plot

unset size; unset ori; unset form;
! rm  plot.gnu.tmp2
set term x11

!epstopdf GLSRV.ps
#!echo "kpdf GLSRV.pdf &"

