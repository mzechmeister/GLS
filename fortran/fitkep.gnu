set fit errorvariables
set macro

# load last fit parameter
load "GLS.log"

# calculate mean weights
wmean=1;
@chi fit wmean file  us 1:(1/$3**2) via wmean
err=":(1)"      # dummy unit errors 
@chi err=":3"

# define Kepler equation
Ex(x)=Ex_sol(x,x,x+e*sin(x)+e**2/2*sin(2*x))
# solve recursive
Ex_sol(x,F,Fn)=(abs((Fn-F)/F)<0.00001)? Fn:Ex_sol(x,Fn,Fn+(x-Fn+e*sin(Fn))/(1-e*cos(Fn)))

# define radial velocity equation
RV(x,RV0,K,e,w,T0,P) = RV0 + K*(e*cos(w*pi/180) +cos(w*pi/180+2*atan(((1+e)/(1-e))**0.5*tan(Ex((x-T0)/P*2*pi)/2))))

# use differential Periastronpassage
# seems that it works better
dt0=0.0001            # init (exact zero gives singular matrix)

fit RV(x,RV0,K,e,w,T0+dt0,PKep) file  us 1:2@err via PKep,dt0,e,w,K,RV0

T0=T0+dt0
dt0=0.

# output
set print "iGLSKep.log"
print ""
print "file ='",file,"'"
#print "chi2= ",FIT_WSSR
print "rms =  ", sqrt(FIT_WSSR/FIT_NDF/wmean)
@chi print "meanerr = ", sqrt(1/wmean)
print ""

print "PKep  = ",PKep," +/- ",PKep_err
print "T0    = ",T0  ," +/- ",dt0_err
print "e     = ",e   ," +/- ",e_err
print "w     = ",w   ," +/- ",w_err
print "K     = ",K   ," +/- ",K_err
print "RV0   = ",RV0 ," +/- ",RV0_err
set print

! cat iGLSKep.log
@calmass system(sprintf("GLS -M %g %g %g %g",mstar,K,PKep,e))

# create a file that emulates RVKepRes.plt
# JD  RVKepRes      RVerr   PhaseKep         RV
set table "iRVKepRes.plt"
phase(x,P)=x/P-ceil(x/P)+1
plot file us 1:($2-RV($1,RV0,K,e,w,T0,PKep))@err:(phase($1-T0,PKep)):2:(RV($1,RV0,K,e,w,T0,PKep)) w xye

! sed -i '2,4d; s/i$//; 1c # JD  RVKepRes      RVerr   PhaseKep       RV     RV_model' iRVKepRes.plt

#set samp 1000
#set term x11 1
#pl RV(x,RV0,K,e,w,dt0,PKep), file us 1:2:3 w e 3
#set term x11 2
#pl file us 1:($2-RV($1,RV0,K,e,w,dt0,PKep)):3 w e 3
#load 'GLSRV.gnu'
set term x11
reset

#pl [0:1] RV(x*PKep+t0,RV0,K,e,w,t0,PKep),"iRVKepRes.plt" us 4:5:3 w e t ""
#pl [0:1] "iRVKepRes.plt" us 4:2:3 w e t ""

