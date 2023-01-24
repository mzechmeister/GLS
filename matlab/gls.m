function [f,pow,fap] = gls_MZ(t, y, f, e, fap_vals)
% [f,pow,fap] = gls(t, y, f, e, fap_vals)
% Calculates The generalized Lomb-Scargle periodogram (Zechmeister 2009).
% If errors are not given, then unweighted means are used.
%
% Input: t - time vecotor of the data
%        y - the data
%        f - frequencies vector to calculate GLS (optional)
%        e - errors of the data (optional)
%        fap_vals - GLS power values to calculate the fap (default = [0.1 0.01 0.001])
%
% Output: f - frequencies.
%         pow - normalized GLS power for each f.
%         fap - selected false-alarm probability GLS values.
%
% Created: 2015, Shay Zucker
% Modified: 2018, Lev Tal-Or & Mathias Zechmeister

% data
N = length(t);
% center time
tmin = min(t);
t = t - tmin;
tbase = max(t);

% default f:
if nargin<3
    ofac = 10;
    hifac = 1;q
    fstep = 1 / tbase / ofac ; % frequency sampling depends on the time span, default for start frequency
    fbeg = fstep;
    fnyq = 0.5 / tbase * N;     % Nyquist frequency
    fend = fnyq * hifac;
    f = fbeg:fstep:fend;
end

% Allocate output array.
%
pow = 0 * f;

% Convert frequencies to angular frequencies.
%
w = 2 * pi * f;

%--------------------------------------------------------------------------
% calculate the FAP lines (analytic aproximation):
%
frange = max(f) - min(f);
M = frange * tbase;
if nargin<5
    fap_vals = [0.1 0.01 0.001];
end
P = 1 - (1 - fap_vals).^(1/M);
fap = 1 - P.^(2/(N-3));

%--------------------------------------------------------------------------
if nargin<4 % if erors are not given use means
%
% Calculate values that are independent of frequency
%
Y = mean(y);
YY_hat = mean(y.^2);
YY = YY_hat - Y^2;

for k = 1:length(f)

    x = w(k) * t;   % phases

    cosx = cos(x);
    sinx = sin(x);
    cos2x = cos(2*x);
    sin2x = sin(2*x);

    omegatau = atan2(mean(sin2x)-2*mean(cosx)*mean(sinx),...
        mean(cos2x)-(mean(cosx)^2-mean(sinx)^2))/2;

    cosx_tau = cos(x-omegatau);
    sinx_tau = sin(x-omegatau);

    C = mean(cosx_tau);
    S = mean(sinx_tau);

    YC_hat = mean(y.*cosx_tau);
    YS_hat = mean(y.*sinx_tau);

    YC = YC_hat - Y*C;
    YS = YS_hat - Y*S;

    CC_hat = mean(cosx_tau.^2);
    SS_hat = mean(sinx_tau.^2);

    CC = CC_hat - C^2;
    SS = SS_hat - S^2;

    pow(k) = 1/YY * (YC^2/CC+YS^2/SS);

end

%--------------------------------------------------------------------------
else % if erors are given use wmeans
%
% Calculate values that are independent of frequency
%
wei = 1./e.^2;
Y = wmean(y,wei);
YY_hat = wmean(y.^2,wei);
YY = YY_hat - Y^2;

for k = 1:length(f)

    x = w(k) * t;

    cosx = cos(x);
    sinx = sin(x);
    cos2x = cos(2*x);
    sin2x = sin(2*x);

    omegatau = atan2(wmean(sin2x,wei)-2*wmean(cosx,wei)*wmean(sinx,wei),...
        wmean(cos2x,wei)-(wmean(cosx,wei)^2-wmean(sinx,wei)^2))/2;

    cosx_tau = cos(x-omegatau);
    sinx_tau = sin(x-omegatau);

    C = wmean(cosx_tau,wei);
    S = wmean(sinx_tau,wei);

    YC_hat = wmean(y.*cosx_tau,wei);
    YS_hat = wmean(y.*sinx_tau,wei);

    YC = YC_hat - Y*C;
    YS = YS_hat - Y*S;

    CC_hat = wmean(cosx_tau.^2,wei);
    SS_hat = wmean(sinx_tau.^2,wei);

    CC = CC_hat - C^2;
    SS = SS_hat - S^2;

    pow(k) = 1/YY * (YC^2/CC+YS^2/SS);

end

end
