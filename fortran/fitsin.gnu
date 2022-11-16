set fit errorvariables
set macro

# load last fit parameter
load "GLS.log"

# calculate mean weights
wmean=1;
@chi fit wmean file  us 1:(1/$3**2) via wmean;
err=":(1)"      # dummy unit errors 
@chi err=":3"

# define radial velocity equation
RV(x,Coff,A,Tph0,P) = Coff + A*sin(2*pi*(x-Tph0)/P)

dTph0=0.0001            # init (exact zero gives singular matrix)

fit RV(x,Coff,Amp,Tph0+dTph0,PSin) file  us 1:2@err via Amp,Coff,PSin,dTph0

Tph0=Tph0+dTph0
dTph0=0.

# output
set print "iGLSSin.log"
print ""
print "file ='",file,"'"
#print "chi2= ",FIT_WSSR
print "rms =  ", sqrt(FIT_WSSR/FIT_NDF/wmean)
@chi print "meanerr = ", sqrt(1/wmean)
print ""

print "PSin  = ",PSin," +/- ",PSin_err
print "Tph0  = ",Tph0," +/- ",dTph0_err
print "Amp   = ",Amp ," +/- ",Amp_err
print "Coff  = ",Coff," +/- ",Coff_err
set print

! cat iGLSSin.log
@calmass system(sprintf("GLS -M %g %g %g %g",mstar,Amp,PSin,0.))

# create a file that emulates RVSinRes.plt
# JD  RVSinRes      RVerr   PhaseSin         RV
set table "iRVSinRes.plt"
phase(x,P)=x/P-ceil(x/P)+1
plot file us 1:($2-RV($1,Coff,Amp,Tph0,PSin))@err:(phase($1-Tph0,PSin)):2:(RV($1,Coff,Amp,Tph0,PSin)) w xye

! sed -i '2,4d; s/i$//; 1c # JD  RVSinRes      RVerr   PhaseSin       RV     RV_model' iRVSinRes.plt

#set samp 1000
#pl file,RV(x,Coff,Amp,Tph0,PSin)
