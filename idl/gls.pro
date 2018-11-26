;+
; Name        : GLS
; Version     : v1.04 (2012/11/08)
;
; Purpose     : Calculate the generalised Lomb-Scargle periodogram
;                                (Zechmeister & Kuerster, 2009, A&A, 496, 577)
;
; Created by  : M. Zechmeister, T. Reinhold
;
; Example     : IDL> t = findgen(100)
;               IDL> y = 30*cos(2.0*!dpi*(t-30)/27.0)+10
;               IDL> gls, f,p, t,y
;
; Syntax      : IDL> gls, f, p, time, ydata [, dy
;                       fbeg=fbeg, pend=pend, fend=fend, pbeg=pbeg,
;                       ofac=ofac, ind=ind,
;                       /ls, /nostat, /noplot, /refine]
;
; Inputs      : t  - time values
;               y  - data values
;               dy - data errors (if given, then weighted periodogram)
;
; Outputs     : f  - frequency array
;               p  - periodogram power array
;
;               [optional outputs]
;               yres   - residuals of data
;               sinpar - parameters of the sine wave
;                        [frequency, amplitude, phase T0, offset]
;                        hint: plot, sinewave(time, sinpar)
;
; Keywords    : fbeg   - start frequency (default: fbeg = 1./timebase/ofac)
;
;               fend   - end frequency (default: fend = .5*N/timebase
;                                        = average Nyquist frequency)
;
;               pbeg   - start period (default: pbeg = timebase*ofac)
;                        ignored, if fend is given
;
;               pend   - end period (default: pend = timebase*ofac)
;                        ignored, if fbeg is given
;
;               ls     - use the classical Lomb-Scargle periodogram
;
;               ofac   - oversampling factor (default: ofac = 20)
;
;               nostat - no statistical info on screen
;
;               noplot - no graphical output
;
;               refine - uses the GLS solution as a starting guess for a
;                        Levenberg-Marquardt fit to obtain the best parameters
;
; Tips:  Prewhitening
;        gls, f,p, t,y    ,yres=yres1
;        gls, f,p, t,yres1,yres=yres2
;        ...
;-
; changes:
; bugfix T0
; output of documentation for invalid nparams()
; code cleaned
; one plot with three panel (instead of the three windows)

FUNCTION sinewave, x,params
;  params=[Freq, Amp, T0, offset]
   nsin = n_elements(params)/4
   par = reform(params,4,nsin)
   y = 0
   FOR i=0,nsin-1 DO y += par[1,i]*SIN(2*!DPI*par[0,i]*(x-par[2,i]))+par[3,i]
   RETURN, y
END

FUNCTION fitsin, t,y,params,dy=dy, yres=yres,BESTNORM=BESTNORM
; improves the fit parameters using MPFIT
; params can be an array of sinpar (see example below) and must be
; an array, dblarr(4) or dblarr(*,4), each row containing [Freq, Amp, T0, offset]
; can fit simultaneously multiple sine waves
; Example     : IDL> t = findgen(100)
;               IDL> y = 30*sin(2.0*!dpi*(t-30.98)/27.0)+10
;               IDL> y = y + 10*sin(2.0*!dpi*(t-17.23)/17.0)-5
;               IDL> gls, f,p, t,y,    sinpar=par1,yres=yres
;               IDL> gls, f,p, t,yres, sinpar=par2
;               IDL> plot, t, y,ps=4 & oplot, t, sinewave(t, [[par1],[par2]]), col=244
;               IDL> parbest = fitsin(t,y,[[par1],[par2]])     ; simultaneous fit
;               IDL> oplot, t, sinewave(t,parbest), col=150
;               IDL> print, ["Freq  Amp T0 offset"],par1,par2, ["new"], parbest

   nsin = n_elements(params)/4
   sinpar_new = params
   sinpar_new[3,*] = 0                      ; set all offsets to zero
   sinpar_new[3,0] = total(params[3,*])     ; combine to a total offset
   parinfo = REPLICATE({fixed:0}, nsin*4)
   IF (nsin GT 1) THEN parinfo[7+4*INDGEN(nsin-1)].fixed = 1.  ; fix all offsets, but not the first
   IF NOT KEYWORD_SET(dy) THEN dy = 1+0.*y  ; unit weights
   sinpar_new = MPFITEXPR('sinewave(x,P)', t,y,dy, sinpar_new,parinfo=parinfo, /QUIET,BESTNORM=BESTNORM)
   IF KEYWORD_SET(yres) THEN yres = y-sinewave(t, sinpar_new)
   RETURN, sinpar_new
END



; ######## main program ######################

PRO GLS, f, p, time, ydata, dy,$
         fbeg=fbeg, fend=fend,  pbeg=pbeg, pend=pend, ofac=ofac, ls=ls,$
         nostat=nostat, noplot=noplot, refine=refine,$
         sinpar=sinpar, yres=yres

IF N_PARAMS() LT 4 THEN BEGIN ; output gls.pro documentation
   OPENR, lun, (routine_info('GLS',/source)).path, /GET_LUN
   glsdoc=''
   WHILE (glsdoc NE ";-") DO BEGIN
      READF, lun, glsdoc
      PRINT, glsdoc
   ENDWHILE
   FREE_LUN, lun
   PRINT, (routine_info('GLS',/source)).path
   RETALL
ENDIF

twopi = 2.0*!DPI
IF NOT KEYWORD_SET(ofac) THEN ofac = 20.    ; set default value

ind = WHERE(FINITE(ydata),nt)      ; check for NaN values, nt = number of valid data points
;                        (fitsin/MPFIT works only with valid data points)
t_data = time[ind]                 ; keep only the valid data points
y_data = ydata[ind]


; prepare input
t = t_data - MIN(t_data)        ; shift the time serie
tbase = MAX(t)                  ; time base

IF N_ELEMENTS(dy) THEN y_err = dy[ind]     ; keep weights of valid data points
w = N_ELEMENTS(y_err)? 1./(y_err)^2 : DOUBLE(1.+0.*y_data)     ; use user or unit weights
w /= TOTAL(w)                ; normalize weights, now TOTAL(w)=1

Y  = TOTAL(w*y_data)        ; Eq.(7), Y is also the (weighted) mean of ydata
wy = y_data-Y               ; shifting data to zero mean
;Y  = TOTAL(w*wy)            ; Y should be now zero
YY = TOTAL(w*wy^2)          ; Eq.(10), variance for the zero mean data
wy = w*wy                   ; data with attached weights

delta_f = 1./tbase/ofac     ; frequency sampling depends on the time span
fmin    = delta_f           ; default for start frequency
delta_t = tbase/(nt-1)      ; mean time sampling
fmax    = 0.5/delta_t       ; mean Nyquist frequency as default for end frequency


; check if the user has specified other ranges
IF NOT KEYWORD_SET(fbeg) THEN fbeg = KEYWORD_SET(pend)? 1./pend : fmin
IF NOT KEYWORD_SET(fend) THEN fend = KEYWORD_SET(pbeg)? 1./pbeg : fmax

; check frequency limits
IF (fend LT fbeg) THEN BEGIN
   PRINT, 'WARNING: fend < fbeg; swapping both values!'
   fmin=fend & fend=fbeg & fbeg=fmin
ENDIF

nout = LONG((fend-fbeg)/delta_f)   ; size of frequency grid

; general infomation
IF NOT KEYWORD_SET(nostat) THEN BEGIN
   PRINT,'  number of input points:     ', nt
   PRINT,'  weighted mean of dataset:   ', Y
   PRINT,'  weighted rms of dataset:    ', SQRT(YY)
   PRINT,'  time base:                  ', tbase
   PRINT,'  number of frequency points: ', nout
ENDIF

f = FINDGEN(nout)*delta_f+fbeg            ; set up the frequency grid
P = DOUBLE(0*f) & A=P & B=P & off=P       ; initialize the power grid

; calculate the periodogram
FOR k=0L,nout-1 DO BEGIN
   omega = twopi*f[k]    ; angular frequency

   cosx = COS(omega*t)   ; to call the cos-function only one time
   sinx = SIN(omega*t)
   C = TOTAL(w*cosx)     ; Eq.(8)
   S = TOTAL(w*sinx)     ; Eq.(9)
   IF KEYWORD_SET(ls) THEN BEGIN ; set C=S=0 for the classical Lomb-Scargle
      C=0 & S=0
   ENDIF
   YC = TOTAL(wy*cosx)   ;   & YC=YC-Y*C     ; Eq.(11) simplifies, since Y should be 0 (mean was subtracted)
   YS = TOTAL(wy*sinx)   ;   & YS=YS-Y*S     ; Eq.(12)   -"-
   CC = TOTAL(w *cosx^2)                     ; Eq.(13)
   SS = 1.-CC                & SS=SS-S*S     ; Eq.(14)
   CC = CC-C*C                               ; Eq.(13)
   CS = TOTAL(w *cosx*sinx)  & CS=CS-C*S     ; Eq.(15)
   D  = CC*SS-CS^2                           ; Eq.(6)

   A[k] = (YC*SS-YS*CS) / D
   B[k] = (YS*CC-YC*CS) / D
   off[k] = -A[k]*C-B[k]*S

   p[k] = (SS*YC^2/D+CC*YS^2/D-2.*CS*YC*YS/D)/YY    ; Eq.(5)
ENDFOR


; preparing output
pmax = MAX(p,i)                           ; returns also the index i of the maximum power
rms = SQRT(YY*(1-p[i]))
; get the curvature in the power peak by fitting a parabola y=aa*x^2
edgeflag = 0
IF (i GT 0 AND i LT nout-1) THEN BEGIN
   xh = (f[i-1:i+1]-f[i])^2               ; shift the parabola origin to power peak
   yh = p[i-1:i+1]-p[i]
   aa = TOTAL(yh*xh)/TOTAL(xh*xh)         ; calculate the curvature (final equation from least square)
   f_err    = ' +/- ' + STRTRIM(SQRT(-2./nt* p[i]/aa*(1-p[i])/p[i]),1)
   Psin_err = ' +/- ' + STRTRIM(SQRT(-2./nt* p[i]/aa*(1-p[i])/p[i])/f[i]^2,1)
ENDIF ELSE BEGIN
   edgeflag = 1
   f_err = ""
   Psin_err = ""
   PRINT, 'WARNING: Highest peak is at the edge of the frequency range!!!'
   PRINT, '   No output of frequency error.'
   PRINT, '   Increase frequency range to sample the peak maximum.'
ENDELSE

fbest = f[i]
Amp = SQRT(A[i]^2+B[i]^2)
ph  = ATAN(A[i],B[i]) / twopi
T0  = MIN(t_data)-ph/f[i]
offset = off[i]+Y                                  ; re-add the mean
sinpar = [fbest,Amp,T0,offset]                     ; summarise

; statistics:
IF NOT KEYWORD_SET(nostat) THEN BEGIN
   PRINT,''
   PRINT,'  maximum power p:        ', p[i]
   PRINT,'  rms of residuals:       ', SQRT(YY*(1-p[i]))
   IF KEYWORD_SET(dy) THEN $
      PRINT,'  mean weighted internal error: ', SQRT(nt/TOTAL(1./y_err^2))
   PRINT,'  best sine frequency f: ', fbest, f_err
   PRINT,'  best sine period Psin: ', 1./fbest, Psin_err
   PRINT,'   amplitude: ', Amp,   ' +/- ',STRTRIM(SQRT(2./nt)*rms,1)
   PRINT,'   phase ph:  ', ph,    ' +/- ',STRTRIM(SQRT(2./nt)*rms/Amp/twopi,1)
   PRINT,'   phase T0:  ', T0,    ' +/- ',STRTRIM(SQRT(2./nt)*rms/Amp/twopi/f[i],1)
   PRINT,'   offset:    ', offset,' +/- ',STRTRIM(SQRT(1./nt)*rms,1)
;   PRINT, '   Psin_err:', SQRT(6./nt)*2/twopi/tbase*rms/Amp/f[i]^2
;   PRINT, '   f_err:   ', SQRT(6./nt)*2/twopi/tbase*rms/Amp
;   PRINT, '   f_err:   ', SQRT(-2./nt* p[i]/aa*(1-p[i])/p[i]) ; error from peak curvature
ENDIF


IF KEYWORD_SET(refine) THEN BEGIN
   sinpar_new = MPFITEXPR('sinewave(x,P)', t_data, y_data, y_err, sinpar, /QUIET,PERROR=perr,BESTNORM=BESTNORM,DOF=DOF)
   perr *= SQRT(BESTNORM / DOF)      ; rescale errors with sqrt(chi^2_red)
   yres = y_data - sinewave(t_data,sinpar_new)
   sinpar = sinpar_new  ; update the returned parameters
   fbest = sinpar[0]
   Amp = sinpar[1]
   T0 = sinpar[2]
   offset = sinpar[3]
   PRINT,'refined parameters:'
   PRINT,'  rms of residuals:      ', SQRT(TOTAL(yres^2*w))
   PRINT,'  best sine frequency f: ', fbest, ' +/- ',STRTRIM(perr[0],1)
   PRINT,'  best sine period Psin: ', 1./fbest, ' +/- ',STRTRIM(perr[0]/fbest^2,1)
   PRINT,'   amplitude: ', Amp,    ' +/- ',STRTRIM(perr[1],1)
   PRINT,'   phase T0:  ', T0,     ' +/- ',STRTRIM(perr[2],1)
   PRINT,'   offset:    ', offset, ' +/- ',STRTRIM(perr[3],1)
ENDIF

; calculate the residuals
yres = ydata - sinewave(time,sinpar)

; graphical output
IF NOT KEYWORD_SET(noplot) THEN BEGIN
   loadct,39
   ; plot periodogram (power(frequency))
   !p.multi=[0,1,3]
   PLOT, f,p, xtitle="frequency", ytitle="power"
   IF NOT edgeflag THEN $
      OPLOT, f,p[i]+aa*(f-f[i])^2,color=244  ; peak curvature

   ; plot time serie (t_i,y_i)
   PLOT, t_data,y_data,psym=4, xtitle="time"
   IF N_ELEMENTS(y_err) THEN OPLOTERR, t_data,y_data, y_err, 4
   t_x = MIN(t_data)+tbase*(FINDGEN(200*tbase*fbest)/(200*tbase*fbest))
   OPLOT, t_x, sinewave(t_x,sinpar), color=244
   t_x = FINDGEN(1000)/1000
   x = Amp*sin(twopi*t_x)+offset

   ; plot phase folded curve
   PLOT, ((t_data-T0) MOD (1./fbest))*fbest,y_data, psym=4, xtitle="phase", xra=[0,1]
   IF N_ELEMENTS(y_err) THEN OPLOTERR, ((t_data-T0) MOD (1./fbest))*fbest,y_data, y_err, 4
   OPLOT, t_x,x,color=244
   WSHOW
   !p.multi=0
ENDIF

END ; of program

