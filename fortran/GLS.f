c     ---------------------------------------------------------------
c     GLS version 2.3.02
c     Calculates the generalised Lomb-Scargle periodogram
c                                  (Zechmeister & Kuerster, 2009, A&A, 496, 577)
c     This program searches for periods with help of
c     the best sine fit and the best Keplerian fit.
c     by Mathias Zechmeister 12.09.2011 zechmeister@astro.physik.uni-goettingen.de
c     2009-05-01: improved e-T0 grid
c     2009-07-29: checking for missing error values; phase [0:1] corrected; typo corrected
c     2010-03-15: parameters output for gnuplot
c                 x0 saved as T0
c                 progress indication
c                 Fix read of empty lines
c     2010-05-21: bug in T0 fixed
c                 output weighted rms
c                 output analytical FAP for Keplerian orbits
c                 optional fast computation for Keplerian periodogram
c     2010-08-14: error estimation
c                 new command line options
c     2010-11-25: unbiased rms for error estimation
c     2011-01-24: correct rms_0 output
c                 ignores empty lines and comment lines containing '#'
c     2011-04-11: input options -k (replaces nk) and -m
c     2011-06-20: code cleaned and formatted
c                 suggests for start and end frequency
c     2011-09-13: output of mass function
c                 option -M
c     2012-08-03: fixed: precision of pi
c                 changed: computation of Amp = sqrt(A^2+B^2)
c     ---------------------------------------------------------------

      SUBROUTINE SineFit(A,B,off,pow,powLS)
c     fit sine wave y=A*cosx+B*sinx+off
c     A,B,off - fit parameter
      DOUBLE PRECISION pow,powLS,A,B,off, CC,SS,CS,C,S,YC,YS,YY,D
      DOUBLE PRECISION v(10000),z(10000),ww(10000),cosx,sinx
      COMMON /data/ v,z,ww,YY,N  ! v = phase, z = data, ww = weights
      CC = 0.
      CS = 0.
      C  = 0.
      S  = 0.
      YC = 0.
      YS = 0.
      do 20 i=1,N                        ! calculate the weighted sums
         cosx = cos(v(i))
         sinx = sin(v(i))
         CC = CC + ww(i) * cosx**2
         CS = CS + ww(i) * cosx*sinx
         C  = C  + ww(i) * cosx
         S  = S  + ww(i) * sinx
         YC = YC + z(i) * cosx
         YS = YS + z(i) * sinx
  20  continue
      SS = 1. - CC
      D  = CC*SS - CS*CS
      powLS = (SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D) / YY         ! Lomb-Scargle power
      CC = CC - C*C
      SS = SS - S*S
      CS = CS - C*S
      D  = CC*SS - CS*CS
      A = (YC*SS-YS*CS) / D
      B = (YS*CC-YC*CS) / D
      off = -A*C-B*S
c      pow= (A*YC+B*YS)/ YY
      pow = (SS*YC**2/D+CC*YS**2/D-2*CS*YC*YS/D) / YY           ! GLS power
      END

      SUBROUTINE Spectral_Window(v,N,WS)
      INTEGER i,N
      DOUBLE PRECISION WS,WC,v(10000)
      WC = 0.
      WS = 0.
      do 30 i=1,N
         WC = WC + cos(v(i))
         WS = WS + sin(v(i))
  30  continue
      WS = (WC**2+WS**2) / N**2
      return
      END

      FUNCTION TrAn(M,e)                   ! calculation of the true anomaly
      DOUBLE PRECISION TrAn,M,Mn,e,F,Fn
      INTEGER i,itermax/8/
      Fn = M+e*sin(M)+e**2/2*sin(2*M)      ! initial value
c     iterative solving of the transcendent Kepler's equation
      do 40 i=1,itermax
         F  = Fn
         Mn = F-e*sin(F)
         Fn = F+(M-Mn)/(1-e*cos(F))
         if (abs((Fn-F)/F) .lt. 0.00001) goto 41
 40   continue
 41   TrAn = 2.0*atan(sqrt((1.+e)/(1.-e))*tan(Fn/2))
      return
      END

      FUNCTION TrAnfast(M,en)              ! fast calculation of the true anomaly
      DOUBLE PRECISION TrAnfast,M,EeM(1:100,0:1024),frac,phase
      INTEGER en,Mn
      COMMON /TrueAnInit/ EeM
      phase = mod(mod(M,1.D0)+1.,1.D0)*1024.
      Mn = aint(phase)
      frac = phase-Mn
c      write(6,*) EeM(en,Mn),EeM(en,Mn+1),(EeM(en,Mn)+EeM(en,Mn+1))*.5
      TrAnfast = EeM(en,Mn) + (EeM(en,Mn+1)-EeM(en,Mn))*frac    ! interpolate
      return
      END

      SUBROUTINE SetupTrAneM(emin,emax,estep)                   ! setup interpolation table
      DOUBLE PRECISION v(1:100,0:1024),TrAn,M,e,twopi,pi
     1  emin,emax,estep
      INTEGER en/0/,nM
      COMMON /CONST/ twopi
      COMMON /TrueAnInit/ v
      pi = DATAN2(0.d0,-1.d0)
c      open(unit=3, file='GLSEinit.tab')
      do 50 e=emin,emax,estep
         en = en+1
         do 51 nM=0, 1024
            M = nM/1024.*twopi
            v(en,nM) = TrAn(M,e)
            if (M .gt. pi) v(en,nM) = v(en,nM)+twopi
            if (v(en,nM) .lt. 0) v(en,nM) = v(en,nM)+twopi
c            write(3,*) M,e,v(en,nM) ,en,nM
  51     continue
  50  continue
c      close(3)
      END

      SUBROUTINE Phases(JD,twopiFreq,N,phase)
      DOUBLE PRECISION JD(*),twopiFreq,phase(*)
      INTEGER i,N
      do 60 i=1,N
         phase(i) = JD(i) * twopiFreq
  60  continue
      END

      SUBROUTINE LIES(text,variab)       ! allows piping of command line arguments
      INTEGER k/1/,iargc
      SAVE k
      CHARACTER text*(*), variab*(*)
      write(6,'(a,$)') text                 ! output of prompt text
      if (iargc() .ge. k) then              ! if in argument list
         call getarg(k,variab)              ! read out here
         write(6,*) variab(:lnblnk(variab)) ! output the value
         k = k+1                            ! position to next argument
      else
         read(5,'(a)') variab               ! else read manually
      endif
      END

      SUBROUTINE LIESd(text,value)       ! LIES for double precision
      DOUBLE PRECISION value
      CHARACTER text*(*), varia*100
      call LIES(text,varia)
c      read(varia,*) value
      if (lnblnk(varia).gt.0) read(varia,*) value  ! default value
      if (lnblnk(varia).eq.0) write(6,*) value
      END

      SUBROUTINE masscalc(mass,K,P,e,err_fac,mass2)
      DOUBLE PRECISION mass,K,P,e, mass2,A1sini,a_axis,
     1   twopi,G,AU,Msun,Mjup,err_fac,massfunc
      PARAMETER(G=6.67384E-11, AU=1.49598E11, Msun=1.98892E30,
     1          Mjup=Msun/1047.56)
      COMMON /CONST/ twopi
      mass2 = 0.
      if (e .eq. 0.) print *," Companion minimum mass Sine fit"
      if (e .ne. 0.) print *," Companion minimum mass Kepler fit"
      massfunc = 86400./twopi/G/Msun*P*(ABS(K)*sqrt(1.-e*e))**3
      print *, " mass function [Msun]:", massfunc
c     iteration to determine the companion mass m2
      do 70 i=1,10
         mass2 = ((mass+mass2)**2*massfunc)**(1./3)
  70  continue
      A1sini = K/twopi*86400.*P*sqrt(1-e**2)
      a_axis = (mass+mass2)/mass2 *A1sini
      print *," Mmin  [Mjup]:  ", mass2*Msun/Mjup
     1           ," +/-", err_fac*mass2*Msun/Mjup/3.
      print *," Mmin  [Msun]:  ", mass2
     1           ," +/-", err_fac*mass2/3.
      print *," a1*sini [AU]:  ", A1sini/AU
     1           ," +/-", err_fac*A1sini/AU
      print *," a       [AU]:  ", a_axis/AU
     1           ," +/-", err_fac*a_axis/AU
      mass2 = mass2*Msun/Mjup
      END

c     ------------------------------------------------------------------
      PROGRAM Main
      implicit none
      INTEGER i,Nmax,N,wexp/2/,io_error/0/,index,iargc,en ! wexp = weightening exponent (default Chi2)
c     !!! weightening of errors: wexp=0 (variance), wexp=2 (Chi2) !!!
      LOGICAL exists, out/.FALSE./,calkep/.FALSE./,calsin/.TRUE./,
     1  fast/.FALSE./,calmass/.FALSE./,onlymass/.FALSE./
      CHARACTER TV*40, name*600, fittype*8/"unknown"/,line*1024
      PARAMETER(Nmax=10000)  ! max. data points
      DOUBLE PRECISION twopi,
     1   JD(Nmax),RV(Nmax),RVerr(Nmax),phase(Nmax),
     2   t(Nmax), wy(Nmax),ww(Nmax), v(Nmax),
     3   JDmin/1.E38/,tbase/0/,RVmean/0/,wsum,YY/0/,
     4   rmssin/0/,rmskep/0/,chisin/0/,chikep/0/,wrmssin/0/,wrmskep/0/,
     5   pow,powLS,powSin/0/,powKep/0/,powKe,SW,
     6   Freq,Fbeg,Fend,step,ofac/20./,
     7   ee, emin/0./, emax/.91/,  estep/.1/,          ! default values
     8   x0, x0min/0./,x0max/359./,x0step,             ! default values
     9   A,B,C, PSin,Amp,ph,CBest, PKep,K,RV0,T0,e,w,
     1   dummy,err_fac, z,TrAn,TrAnfast,dRVsin,dRVkep,
     2   mass,mass2S,mass2K,M,FAP,prob,FAPKep,probKep
      COMMON /data/ v,wy,ww,YY,N
      COMMON /CONST/ twopi
      twopi = 2.*DATAN2(0.d0,-1.d0)

c     read parameters
      inquire(file="GLS.par", exist=exists)
      if (exists) then
         write(6,*) "found parameter file GLS.par!!!"
         open(unit=2, file= "GLS.par", status= "OLD")
         read(2,*) ofac                         ! default oversampling factor = 20
         read(2,*) wexp                         ! default wexp = 0 (variance)
         read(2,*) emin,emax,estep              ! default 0.  0.91  0.1
         read(2,*) x0min,x0max                  ! default 0.  359.  in degree
         close(2)
      endif                                     ! (else initialize the defaults)
      x0min = x0min /360.*twopi                 ! convert to rad
      x0max = x0max /360.*twopi

c     checking command line for options -kKvcfmM
      if (iargc() .gt. 0) then           ! if arguments
         i=0
  100    i=i+1
         call getarg(i,name)             ! testing i-th command line argument
         if (index(name,'-') .eq. 1) then
            call LIES(' found option:  ', name)         ! get the argument
            if (index(name,'k') .gt. 0) calkep=.TRUE.   ! i.e. calculate Keplerian periodogram
c            if(index(name,'K').gt.0) calsin=.FALSE.    ! i.e. only Keplerian periodogram (no GLS)
            if (index(name,'K') .gt. 0) calkep=.TRUE.   ! -K needed by user?
            if (index(name,'v') .gt. 0) wexp=0          ! i.e. variance
            if (index(name,'c') .gt. 0) wexp=2          ! i.e. chi2
            if (index(name,'f') .gt. 0) fast=.TRUE.     ! i.e. fast Keplerian computation
            if (index(name,'m') .gt. 0) calmass=.TRUE.  ! i.e. calculate companion mass
            if (index(name,'M') .gt. 0) onlymass=.TRUE. ! i.e. calculate only companion mass
            GOTO 100
         endif
      endif

      if (onlymass) then
         call LIESd('\nmass     [Msun]:  ',mass)
         call LIESd(  'amplitude [m/s]:  ',K)
         call LIESd(  'period      [d]:  ',PKep)
         call LIESd(  'eccentricity   :  ',e)
         call masscalc(mass,K,PKep,e,0.d0,mass2K)
         GOTO 500      ! exit program
      endif

      if (fast) call SetupTrAneM(emin,emax,estep)
      if (fast) write(6,*) "setup done"

      call LIES('\nfilename:  ', name)

c     read data from file and calculate some sums
  110 open(unit=1, file=name, status='old' )
      if (index(name,'.vels ') .gt. 0)  read(1,*)     ! skip header
      if (wexp .eq. 2) fittype = "chi2"
      if (wexp .eq. 0) fittype = "variance"
      wsum = 0.
      do 130 i=1,Nmax
  120    read(1,'(a)',END=140) line                   ! read line from datafile
         if (index(line,'#') .gt. 0) GOTO 120         ! skip comment line
         if (lnblnk(line)    .eq. 0) GOTO 120         ! skip empty line
         if (index(name,'.vels ') .gt. 0) then
            read(line,*) TV,JD(i),RV(i),dummy,RVerr(i)
         elseif (wexp .gt. 0) then                    ! with error values, chi2
            read(line,*,IOSTAT=io_error) JD(i),RV(i),RVerr(i)
c            write(6,*) i,JD(i),RV(i),RVerr(i)
            if (io_error .ne. 0 .or. RVerr(i) .eq. 0.) THEN
               write(6,*) "Forcing to wexp=0 and reopen data file!!!"
               if (RVerr(i).eq. 0.) io_error=-2
               close(1)
               wexp = 0
               GOTO 110  !  to reopen file
            endif
         else                              ! no errors values, variance
            read(line,*) JD(i),RV(i)
            RVerr(i) = 1.
         endif
         JDmin = MIN(JDmin,JD(i))
         N = i
         ww(i) = (1./RVerr(i))**wexp
         wsum  = wsum + ww(i)
 130  continue
 140  close(1)

      do 150 i=1,N
         t(i) = JD(i) - JDmin              ! shift the time series
         tbase = MAX(tbase,t(i))           ! time base
         ww(i) = ww(i) / wsum              ! normalize weights
         RVmean = RVmean + RV(i)*ww(i)     ! weighted mean
 150  continue

      do 160 i=1,N
         wy(i) = RV(i)-RVmean              ! centering RVs
         YY    = YY + wy(i)**2*ww(i)       ! sum for chi2 above the mean
         wy(i) = wy(i)*ww(i)               ! attach weights to centered data (saves many multiplications)
 160  continue

      print *, ".........", N, " data points readed"
      print *, "......... time base ", tbase, " days"
      print *, "......... 1/time base ", 1./tbase, " 1/d"
      print *, "......... ",fittype(:lnblnk(fittype))//"_0", YY*wsum
      print *, "......... rms_0 ",sqrt(YY*N/(N-1))
c     read search intervals
      Fbeg = 1./tbase/10.              ! default value
      Fend = N*Fbeg*10.                ! default value (average Nyquist frequency)
      write(TV,'(G12.5)') Fbeg
      call LIESd(" start frequency [1/d] ("//TV(1:12)//"): ",Fbeg)
      write(TV,'(G12.5)') Fend
      call LIESd(" end frequency   [1/d] ("//TV(1:12)//"): ",Fend)
      step = 1./tbase/ofac
      write(6,*)"calculate ",int((Fend-Fbeg)/step+1)," frequency steps"

      print *, "\n *** GLS periodogram          ********************"
      print *,    "*** (analysis with Sine Fit) ********************"
      open(unit=10, file='GLS.plt')
      write(10,'(4A18)')  "#frequency","pow","window","powLS"
      do 170 Freq=Fbeg,Fend, step
         call Phases(t,twopi*Freq,N,v)
         call SineFit(A,B,C,pow, powLS)
         call Spectral_Window(v,N,SW)
         write(10,*)  Freq, pow, SW, powLS
         if (pow .GT. powSin) then
            powSin = pow
            ph = MOD(DATAN2(A,B)+twopi, twopi)
            Amp = sqrt(A**2+B**2) ! = A/sin(ph)
            CBest = RVmean+C
            PSin = 1./Freq
            print *, "Period: ", PSin, "\tpow: ", pow           ! output of the improved periods
         endif
 170  continue
      close(10)

      if (calkep) then
      print *, "\n\n *** Keplerian periodogram            ************"
      print *,      "*** (analysis with Keplerian orbits) ************"
c     fit Keplerian orbits RV=A*cosv+B*sinv+C
      open(unit=20, file='GLSKep.plt')
      write(20,'(A14,A11)')  "#frequency", "pow"
      do 180 Freq=Fbeg,Fend, step
         call Phases(t,twopi*Freq,N,phase)
         powKe = 0.
         en = 0
         do 190 ee=emin,emax,estep
            en = en+1
            x0step = twopi/int(twopi*ee/estep+1)               ! polar grid e-x0
            do 200 x0=x0min,x0max,x0step
c              with known x0 and e parameter A, B und C can be minimised,
c              because the true anomaly is now calculable:
               if (fast) then
                  do 205 i=1,N
                     v(i) = TrAnfast((phase(i)-x0)/twopi,en)
 205              continue
               else
                  do 206 i=1,N
                     v(i) = TrAn(phase(i)-x0,ee)
 206              continue
               endif
               call SineFit(A,B,C,pow,powLS)
               powKe = MAX(pow, powKe)
               if (pow .GT. powKep) then
                  powKep = pow
                  PKep = 1./Freq
                  T0 = x0*PKep/twopi
                  e = ee
                  w = MOD(-DATAN2(B,A)+twopi, twopi)
                  K = sqrt(A**2+B**2) ! = -B/sin(w)
                  RV0 = RVmean+C-A*e
                  out = .true.
               endif
 200        continue
 190     continue
         if (out) print *, "Period: ", PKep, "\tpow: ", powKep         ! output of the improved periods
         out = .false.
         write(20,*) Freq, powKe
         write(6,'(f10.2,a,$)') (Freq-Fbeg)/(Fend-Fbeg)*100.," %\r"          ! progress indication
 180  continue
      close(20)
      write(6,*) "    100.00 % done!\n"
      endif  ! calkep

c     writing output files
      open(unit=12, file='RVSinRes.plt')
      if (calkep) open(unit=13, file='RVKepRes.plt')
      open(unit=14, file='GLSFit.plt')
      write(12,'(5(A11))') "#JD","RVSinRes","RVerr","PhaseSin","RV"
      if (calkep)
     1    write(13,'(5(A11))') "#JD","RVKepRes","RVerr","PhaseKep","RV"
      write(14,'(A,2(A10))') "#Phase","RVSin","RVKep"
      call Phases(t,1./PSin,N,v)
      call Phases(t,1./PKep,N,phase)
      do 300 i=1,N                      ! RVSinRes.plt and RVKepRes.plt
         dRVsin = RV(i)-(Amp*sin(t(i)*twopi/PSin+ph)+CBest)
         dRVkep = RV(i)-
     1        (RV0+K*(e*cos(w)+cos(TrAn((t(i)-T0)*twopi/PKep,e)+w)))
         write(12,*) JD(i),dRVsin,RVerr(i)
     1           ,MOD(MOD(v(i),1.D0)+1.,1.D0),RV(i)
         if (calkep) write(13,*) JD(i),dRVKep,RVerr(i)
     1           ,MOD(MOD(phase(i),1.D0)+1.,1.D0), RV(i)
         rmssin = rmssin + dRVsin**2
         rmskep = rmskep + dRVkep**2
         chisin = chisin + (dRVsin/RVerr(i))**2
         chikep = chikep + (dRVkep/RVerr(i))**2
 300  continue
      rmssin = sqrt(rmssin/(N-4))           ! unbiased rms
      rmskep = sqrt(rmskep/(N-6))
      wrmssin = sqrt(chisin/wsum*N/(N-4))
      wrmskep = sqrt(chikep/wsum*N/(N-6))

      do 310 z=0,twopi, twopi/100.           ! GLSFit.plt
         write(14,*) z/twopi, Amp*sin(z+ph)+CBest,
     1               RV0+K*(e*cos(w)+cos(TrAn(z-T0*twopi/PKep,e)+w))
 310  continue
      close(12)
      if (calkep) close(13)
      close(14)

      prob= (1.-powSin)**((N-3.)/2)        ! spectral false alarm probability
      M   = tbase * abs(Fend-Fbeg)         ! number of independent frequencies
      FAP = M * prob
      if (FAP .gt. 0.01) FAP = 1. - (1.-prob)**M

      write (6,*) 'used parameters',ofac, wexp, emin,emax,estep
      print *, " ",fittype(:lnblnk(fittype))//" for constant:",YY*wsum
      print *, " rms for constant: ",sqrt(YY*N/(N-1))
      if (wexp .eq. 2)
     1   print *, " mean weighted internal error:", sqrt(N/wsum)
c     Output of best sine fit RV=Amp*sin(2pi((t-t0)/PSin+ph))+CBest
      print *, "\n Best sine fit values for *********** ",
     1   name(:lnblnk(name))
      print *, " best sine period [d]:       ", PSin
     1   ," +/-",sqrt(6./N)*2/twopi/tbase*wrmssin/Amp*PSin**2
      print *, " amplitude Amp:              ", Amp
     1   ," +/-",sqrt(2./N)*wrmssin
      print *, " phase ph [0..1] to JDmin:   ", ph/twopi
     1   ," +/-",sqrt(2./N)*wrmssin/Amp/twopi
      print *, " T0 [JD] to JDmin:\t      ", JDmin-ph/twopi*PSin+PSin
     1   ," +/-",sqrt(2./N)*wrmssin/Amp/twopi*PSin
      print *, " offset RV0:                 ", CBest
     1   ," +/-",sqrt(1./N)*wrmssin
      write(6,*)
      print *, fittype(:lnblnk(fittype))," of residuals:",
     1   (1.-powSin)*YY*wsum,chisin ! variance after fit
      print *, " weighted rms of residuals:  ", wrmssin
      print *, " rms of residuals:           ", rmssin
      print *, " highest power p[0..1]:      ", powSin
      print *, "               P[0..(N-1)/2]:", powSin*(N-1)/2
      print *, " FAP(spectral)=(1-p)^((N-3)/2):", prob
      print *, " FAP:                        ", FAP

      if (calkep) then
c      probKep = (1+ (N-3)/2 * powKep/(1-powKep)) * (1-powKep)**((N-3)/2)
      probKep = (1+ (N-5.)/2 * powKep) * (1-powKep)**((N-5.)/2)
      FAPKep = M * probKep
      if (FAPKep .gt. 0.01) FAPKep = 1. - (1.-probKep)**M
c     Output of best Keplerian fit
      print *, "\n Best Keplerian fit value for ********** ",
     1   name(:lnblnk(name))
      print *, " best Keplerian period [d]:   ", PKep
     1   ," +/-",sqrt(6./N)*2./twopi/tbase*wrmskep/K*PKep**2
      print *, " RV semi-amplitude K:         ", K
     1   ," +/-",sqrt(2./N)*wrmskep
      print *, " eccentricity e [0..1]        ", e
      print *, " longitude of periastron w [deg]:", w*360/twopi
      print *, " periastron passage T0:       ", T0+JDmin
     1   ," +/-",sqrt(2./N)*wrmskep/K/twopi*PKep
      print *, " offset RV0:                  ", RV0
     1   , " +/-",sqrt(1./N)*wrmskep

      write(6,*)
      print *, fittype(:lnblnk(fittype))," of residuals:",
     1   (1.-powKep)*YY*wsum,chikep ! variance after fit
      print *, " weighted rms of residuals:    ", wrmskep
      print *, " rms of residuals:             ", rmskep
      print *, " highest Keplerian power p:    ", powKep
      print *, " FAP(spectral, see Cumming 04):", probKep
      print *, " FAP:                          ", FAPKep
      endif !calkep

      if (calmass) then
         call LIESd("\nstar mass [sun mass]: ",mass)
         err_fac = sqrt(2./N)*wrmssin/Amp
     1             * sqrt(1.+(sqrt(3.)*2./twopi*PSin/tbase)**2)
         call masscalc(mass,Amp,PSin,0.d0,err_fac,mass2S)
         if (calkep) then
         err_fac = sqrt(2./N)*wrmskep/K
     1             * sqrt(1.+(sqrt(3.)*2./twopi*PKep/tbase)**2)
         call masscalc(mass,K,PKep,e,err_fac,mass2K)
         endif
      endif !calmass


c     add results to GLS.log
      open(unit=3,  file='GLS.log',access='APPEND')
      write(3,*)"\nfile='",name(:lnblnk(name)),"'"
      write(3,*)"PSin=",PSin,";  Amp=",Amp,";"
      write(3,*)" Tph0=",JDmin-ph/twopi*PSin+PSin,"; Coff=",CBest
      if (calmass) write(3,*) " m=",mass2S,"; mstar=",mass
      if (calkep) write(3,*) "PKep=",PKep, "; T0=",T0+JDmin,"; e=",e,";"
      if (calkep) write(3,*) " w=",w*360/twopi,"; K=",K,"; RV0=",RV0
      if (calkep .AND. calmass) write(3,*) " mkep=",mass2K
      close(3)

      write(6,*) "***",N," data points of ",name(:lnblnk(name))
     1           ," processed ***"

      if (io_error .ne. 0)  write(6,*)
     1      "\nWARNING: The fit type was forced to variance (wexp=0)!",
     2      "Check your data file for missing or zero error values,",
     3      "if you like to do chi^2 fitting!!!"
 500  STOP
      END
