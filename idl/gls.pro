;+
; NAME:
;       GLS
;
; VERSION:
;       v2.01 (2016-08-02)
;
; PURPOSE:
;       Calculate the generalised Lomb-Scargle periodogram
;               (Zechmeister & Kuerster, 2009, A&A, 496, 577)
;
; AUTHOR:
;       M. Zechmeister, T. Reinhold
;
; CALLING SEQUENCE:
;       p = gls(t, y, [dy,
;                fbeg=, fend=, pend=, pbeg=, ofac=, hifac=
;                /fast, /ls, /plot, /refine, /verbose,
;                freqs=, ymod=, sinpar=, info=])
;
; INPUTS:
;       t  - time values
;       y  - data values
;       dy - data errors (if given, then weighted periodogram)
;
; OUTPUTS:
;       p  - periodogram power array
;
; OPTIONAL OUTPUTS:
;       freqs  - frequency array
;       ymod   - modelled data values
;       sinpar - parameters of the sine wave
;                [frequency, amplitude, phase T0, offset]
;                hint: plot, sinmod(time, sinpar)
;       info   - structure with statistical information
;                {pmax - maximum power
;                 fopt - optimal frequency
;                 Popt - optimal period
;                 rms  - rms of best fit
;                 rms0 - rms before fit
;                 interr - internal data error}
;
; KEYWORDS:
;       fbeg   - start frequency (default: fbeg = 1/timebase/ofac)
;       fend   - end frequency (default: fend = 0.5*N/timebase*hifac
;                                             = average Nyquist frequency)
;       pbeg   - start period (default: pbeg = 2/N*timebase/hifac)
;                ignored, if fend is given
;       pend   - end period (default: pend = timebase*ofac)
;                ignored, if fbeg is given
;       hifac  - maximum frequency above Nyquist (default: hifac = 1)
;       ofac   - oversampling factor (default: ofac = 10)
;       /fast  - uses trigonometric recurrences
;       /ls    - use the classical Lomb-Scargle periodogram
;       /plot  - graphics for periodogram, time series, residuals, phase folding
;       /refine - uses the GLS solution as a starting guess for a
;                Levenberg-Marquardt fit to obtain the best parameters
;       /verbose - prints statistical info
;
; EXAMPLES:
;       GLS without weights
;
;       IDL> t = findgen(100)
;       IDL> y = 30*cos(2.0*!dpi*(t-30)/27.0)+10
;       IDL> plot, f, gls(t, y, freqs=f)
;
;       Pre-whitening with GLS
;
;       IDL> t = findgen(100)
;       IDL> y = 30*sin(2.0*!dpi*(t-30.98)/27.0)+10
;       IDL> y += 10*sin(2.0*!dpi*(t-17.23)/17.0)-5
;       IDL> p = gls(t, y, sinpar=par1, ymod=yc1, /plot, /verb)
;       IDL> p = gls(t, y-yc1, sinpar=par2, /plot)
;       IDL> plot, t, y,ps=4 & oplot, t, sinmod(t, [[par1],[par2]]), col=244
;       IDL> parbest = sinfit(t, y, [[par1],[par2]])     ; simultaneous fit
;       IDL> oplot, t, sinmod(t, parbest), col=150
;       IDL> print, ["Freq  Amp T0 offset"], par1,par2, ["new"], parbest
;
; MODIFICATION HISTORY:
;       v1.05 - code cleaned
;       v2.00 - sums replaced by dot product
;             - gls now as function
;             - subfunctions renamed
;             - glsplot and glsprint included
;       v2.01 - new: fast option
;             - nostat option replaced by verbose
;-

function sinmod, x, par
; INPUT: par = [Freq, Amp, T0, offset]
   y = total(par[3,*])
   for i=0,n_elements(par)-4,4 do $
      y += par[i+1] * sin(2*!dpi*par[i]*(x-par[i+2]))
   return, y
end

function sinfit, t, y, par, dy=dy, ymod=ymod, chisq=chisq
; improves the fit parameters using MPFIT
; can fit simultaneously multiple sine waves
; params can be an array of sinpar (see example below) and must be
; an array, dblarr(4) or dblarr(*,4), each row containing [Freq, Amp, T0, offset]

   nsin = n_elements(par) / 4
   sinpar_new = par
   if nsin gt 1 then begin
      sinpar_new[3,*] = 0                     ; set all offsets to zero
      sinpar_new[3,0] = total(par[3,*])       ; combine to a total offset
      parinfo = replicate({fixed:0}, 4, nsin)
      parinfo[3,1:*].fixed = 1                ; fix all offsets, but not the first
   endif

   sinpar_new = mpfitexpr('sinmod(x, p)', t, y, dy, sinpar_new, parinfo=parinfo, /quiet, bestnorm=chisq, yfit=ymod)

   return, sinpar_new
end

pro glsprint, par, perr, info
   print, ''
   print, '  maximum power p:        ', info.pmax
   print, '  rms of residuals:       ', info.rms
   if finite(info.interr) then $
      print, '  mean weighted internal error: ', info.interr
   fmt = "(A-19,D0.0,' +/- ',D0.0)"
   print, '   period Psin:', 1./par[0], perr[0]/par[0]^2, form=fmt
   print, '   frequency f:', par[0],perr[0], form=fmt
   print, '   amplitude:', par[1], perr[1], form=fmt
   print, '   phase T0:', par[2], perr[2], form=fmt
   print, '   offset:', par[3], perr[3], form=fmt
end

pro glsplot, t, y, yerr, f, p, sinpar
   fopt = sinpar[0]
   T0 = sinpar[2]

   tvlct, 255, 0, 0, 100   ; red
   !p.multi = [0, 1, 3]

   plot, title='Periodogram (fopt='+strtrim(fopt,2)+', P='+strtrim(1/fopt,2)+')', $
        f, p, xtitle='frequency', ytitle='power', charsiz=2
   oplot, fopt*[1,1], [0,1], color=100   ; mark best freq

   !p.multi = [4, 2, 3]
   ; plot time serie (t_i, y_i)
   plot, t, y, psym=4, charsiz=2, xmargin=[10,-8], ymargin=[0,0], xtickformat='(A1)', ytitle='data'
   if n_elements(yerr) gt 0 then oploterr, t, y, yerr, 4
   tmin = min(t, max=tmax)
   tt = tmin + findgen(200*((tmax-tmin)*fopt>1))/(199*fopt)   ; at least one cycle
   yy = sinmod(tt, sinpar)
   oplot, tt, yy, color=100

   ; plot phase folded curve
   x = ((t-T0)*fopt+1) mod 1.
   plot, x, y, psym=4, xra=[0,1], ymargin=[0,0], charsiz=2, xtickformat='(A1)', ytickformat='(A1)'
   if n_elements(yerr) gt 0 then oploterr, x, y, yerr, 4
   xx = ((tt[0:200-1]-T0)*fopt+1) mod 1.   ; take first cycle
   ind = sort(xx)
   oplot, xx[ind], yy[ind], color=100

   ; residuals time serie
   yres = y - sinmod(t, sinpar)
   plot, t, yres, psym=4, xtitle='time', ytitle='residuals', charsiz=2, xmargin=[10,-8]
   if n_elements(yerr) gt 0 then oploterr, t, yres, yerr, 4
   oplot, [tmin,tmax], [0,0], color=100

   plot, x, yres, psym=4, xra=[0,1], xtitle='phase', charsiz=2, ytickformat='(A1)'
   if n_elements(yerr) gt 0 then oploterr, x, yres, yerr, 4
   oplot, [0,1], [0,0], color=100

   !p.multi = 0
end


; ######## MAIN PROGRAM ######################

function GLS, time, ydata, dy, $
         fbeg=fbeg, fend=fend, pbeg=pbeg, pend=pend, ofac=ofac, hifac=hifac, $
         fast=fast, ls=ls, plot=plot, verbose=verbose, refine=refine, $
         freqs=f, sinpar=sinpar, ymod=ymod, info=info

   on_error, 2

   if n_params() lt 2 then begin ; output gls.pro documentation
      doc_library, 'gls'
      return, 0
   endif

   twopi = 2.0 * !dpi

   fast = keyword_set(fast)
   ls = keyword_set(ls)
   if ~n_elements(ofac) then ofac = 10.
   if ~n_elements(hifac) then hifac = 1.

   ; keep only the valid data points and flatten (for dot-product)
   ind = where(finite(ydata), nt)   ; check for NaN values, nt = number of valid data points
                                    ; (sinfit/MPFIT works only with valid data points)
   t_data = double(time[ind])
   y_data = double(ydata[ind])

   ; prepare input
   tmin = min(t_data)
   t = t_data - tmin           ; shift the time serie
   tbase = max(t)              ; time base

   if n_elements(dy) gt 0 then y_err = dy[ind]   ; keep weights of valid data points
   w = n_elements(y_err) gt 0? 1.d/y_err^2 : 0.d*y_data+1   ; use user or unit weights
   wsum = total(w)
   w = 1/wsum # w              ; normalize weights (now total(w)=1) and prepare dot product

   Y = (w # y_data)[0]         ; Eq.(7), Y is also the (weighted) mean of ydata
   wy = y_data - Y             ; shifting data to zero mean
   ;Y  = total(w*wy)           ; Y should be now zero
   YY = (w # wy^2)[0]          ; Eq.(10), variance for the zero mean data
   wy = w * wy                 ; data with attached weights

   fstep = 1 / tbase / ofac  ; frequency sampling depends on the time span, default for start frequency
   fnyq = 0.5 * (nt-1) / tbase    ; mean time sampling

   ; check if the user has specified other ranges
   if not keyword_set(fbeg) then fbeg = keyword_set(pend)? 1./pend : fstep
   if not keyword_set(fend) then fend = keyword_set(pbeg)? 1./pbeg : fnyq*hifac ; default mean Nyquist frequency for end frequency

   ; check frequency limits
   if fend lt fbeg then begin
      print, 'WARNING: fend < fbeg; swapping both values!'
      fmin = fend & fend = fbeg & fbeg = fmin
   endif

   nk = long((fend-fbeg) / fstep)   ; size of frequency grid

   ; general information
   if keyword_set(verbose) then begin
      print, '  number of input points:     ', nt
      print, '  weighted mean of dataset:   ', Y
      print, '  weighted rms of dataset:    ', sqrt(YY)
      print, '  time base:                  ', tbase
      print, '  number of frequency points: ', nk
   endif

   f = fbeg + fstep*dindgen(nk)           ; set up the frequency grid
   p = dblarr(nk) & A=p & B=p & off=p       ; initialize the power grid

   if fast then $
      ; prepare trigonometric recurrences
      eid = exp(dcomplex(0, twopi*fstep*t)) ; cos(dx)+i sin(dx)


   ; calculate the periodogram
   for k=0L,nk-1 do begin
      if fast then begin
         ; refresh recurrences to avoid error propagation
         if k mod 1000L eq 0 then eix = exp(dcomplex(0, twopi*f[k]*t))
         cosx = double(eix)
         sinx = imaginary(eix)
         eix *= eid             ; increase freq for next loop
      endif else begin
         x = twopi * f[k] * t   ; angular phases
         cosx = cos(x)          ; only one call to costly cos-function
         sinx = sin(x)
      endelse

      C = (w # cosx)[0]                         ; Eq.(8)
      S = (w # sinx)[0]                         ; Eq.(9)
      YC = (wy # cosx)[0]    ;  & YC -= Y * C   ; Eq.(11) simplifies, since Y should be 0 (mean was subtracted)
      YS = (wy # sinx)[0]    ;  & YS -= Y * S   ; Eq.(12)   -"-
      CC = (w # (cosx*cosx))[0]                 ; Eq.(13)
      CS = (w # (cosx*sinx))[0]                 ; Eq.(15)
      SS = 1. - CC                              ; Eq.(14)

      if ~ls then begin      ; else C=S=0 for the classical Lomb-Scargle
         CC -= C * C         ; Eq.(13)
         SS -= S * S         ; Eq.(14)
         CS -= C * S         ; Eq.(15)
      endif
      D = CC*SS - CS^2       ; Eq.(6)

      A[k] = (YC*SS-YS*CS) / D
      B[k] = (YS*CC-YC*CS) / D
      off[k] = -A[k]*C - B[k]*S
      p[k] = (SS*YC^2+CC*YS^2-2.*CS*YC*YS) / D / YY   ; Eq.(5)
   endfor

   ; preparing output
   pmax = max(p, k)                        ; returns also the index k of the maximum power
   rms = sqrt(YY * (1-p[k]))

   ; get the curvature in the power peak by fitting a parabola y=aa*x^2
   edgeflag = k le 0 or k ge nk-1
   if ~edgeflag then begin
      xh = (f[k-1:k+1] - f[k])^2           ; shift the parabola origin to power peak
      yh = p[k-1:k+1] - p[k]
      aa = total(yh*xh) / total(xh*xh)     ; calculate the curvature (final equation from least square)
      ; error estimate from parabolic Taylor expansion around chisq mininum and the criterion:
      ; chisq_red = delta_chisq = 0.5 * chisq'' * (f+f_err - f)^2
      ; f_err^2 = 2 chisq_red / chisq'' = 2 / N * chisq / chisq''
      ;    chisq = chisq0 * (1-p) =>  chisq'' = -chisq0 * p''
      ;    aa = p'' (curvature in p)
      ; f_err^2 = -2 / N  * (1-p) / p''
      f_err = sqrt(-2. / nt * (1-p[k]) / aa)
   endif else begin
      f_err = !Values.F_NAN
      print, 'WARNING: Highest peak is at the edge of the frequency range!!!'
      print, '   No output of frequency error.'
      print, '   Increase frequency range to sample the peak maximum.'
   endelse

   fopt = f[k]
   Amp = sqrt(A[k]^2 + B[k]^2)
   ph = atan(A[k], B[k]) / twopi
   T0 = tmin - ph/f[k]
   offset = off[k] + Y                ; re-add the mean
   sinpar = [fopt, Amp, T0, offset]   ; summarise
   sinperr = [f_err, sqrt(2./nt)*rms, $
                     sqrt(2./nt)*rms/Amp/twopi/f[k], $
                     sqrt(1./nt)*rms]

   info = {k:k, fopt:f[k], Popt:1/f[k], pmax:p[k], rms:rms, interr:!Values.D_NAN}
   if keyword_set(dy) then info.interr = sqrt(nt/wsum)

   ; statistics
   if keyword_set(verbose) then $
      glsprint, sinpar, sinperr, info

   if keyword_set(refine) then begin
      sinpar = MPFITEXPR('sinmod(x,P)', t_data, y_data, y_err, sinpar, /quiet, perror=perr, bestnorm=chisq, dof=DOF)
      perr *= sqrt(chisq / DOF)   ; rescale errors with sqrt(chi^2_red)
      info.rms = sqrt(chisq / wsum)   ; sqrt(total(yres^2*w))
      info.fopt = sinpar[0]
      info.Popt = 1 / info.fopt
      info.pmax = (YY-info.rms^2) / YY
      print, 'refined parameters:'
      glsprint, sinpar, sinperr, info
   endif

   ; graphical output
   if keyword_set(plot) then $
      glsplot, time, ydata, dy, f, p, sinpar

   ; calculate the model
   ymod = sinmod(time, sinpar)

   return, p
end   ; of program
