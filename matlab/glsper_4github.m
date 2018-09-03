function [gls,fap] = glsper_4github(f,t,y,e,fap_vals)
% [gls,fap] = glsper_4github(f,t,y,e,fap_vals)
% Calculates The generalized Lomb-Scargle periodogram (Zechmeister 2009).
% If errors are not given, then unweighted means are used.
% 
% Input: f - frequencies vector to calculate GLS
%        t - time vecotor of the data
%        y - the data 
%        e - errors of the data (optional)
%        fap_vals - GLS power values to calculate the fap (default = [0.1 0.01 0.001])
%
% Output: gls - normalized GLS power for each f.
%         fap - selected false-alarm probability GLS values.
%
% Created: 2015, Shay Zucker
% Modified: 2018, Lev Tal-Or


% Allocate output array.
%
gls = f*0;

% Convert frequencies to angular frequencies.
%
w = 2*pi*f;

%--------------------------------------------------------------------------
% calculate the FAP lines (analytic aproximation):
% 
N = length(t);
trange = max(t)-min(t);
frange = max(f) - min(f);
M = frange*trange;
if nargin<5
    fap_vals = [0.1 0.01 0.001];
end
P = 1 - (1 - fap_vals).^(1/M);
fap = 1-P.^(2/(N-3));

%--------------------------------------------------------------------------
if nargin<4 % if erors are not given use means
%
% Calculate values that are independent of frequency
%
Y = mean(y);
YY_hat = mean(y.^2);
YY = YY_hat - Y^2;

for ifreq = 1:length(f)
    
    this_w = w(ifreq);
    
    coswt = cos(this_w*t);
    sinwt = sin(this_w*t);
    cos2wt = cos(2*this_w*t);
    sin2wt = sin(2*this_w*t);
    
    omegatau = atan2(mean(sin2wt)-2*mean(coswt)*mean(sinwt),...
        mean(cos2wt)-(mean(coswt)^2-mean(sinwt)^2))/2;
    
    coswt_tau = cos(this_w*t-omegatau);
    sinwt_tau = sin(this_w*t-omegatau);
        
    C = mean(coswt_tau);
    S = mean(sinwt_tau);
   
    YC_hat = mean(y.*coswt_tau);
    YS_hat = mean(y.*sinwt_tau);
    
    YC = YC_hat - Y*C;
    YS = YS_hat - Y*S;
    
    CC_hat = mean(coswt_tau.^2);
    SS_hat = mean(sinwt_tau.^2);
    
    CC = CC_hat - C^2;
    SS = SS_hat - S^2;
    
    gls(ifreq) = 1/YY*(YC^2/CC+YS^2/SS);
    
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

for ifreq = 1:length(f)
    
    this_w = w(ifreq);
    
    coswt = cos(this_w*t);
    sinwt = sin(this_w*t);
    cos2wt = cos(2*this_w*t);
    sin2wt = sin(2*this_w*t);
    
    omegatau = atan2(wmean(sin2wt,wei)-2*wmean(coswt,wei)*wmean(sinwt,wei),...
        wmean(cos2wt,wei)-(wmean(coswt,wei)^2-wmean(sinwt,wei)^2))/2;
    
    coswt_tau = cos(this_w*t-omegatau);
    sinwt_tau = sin(this_w*t-omegatau);
        
    C = wmean(coswt_tau,wei);
    S = wmean(sinwt_tau,wei);
   
    YC_hat = wmean(y.*coswt_tau,wei);
    YS_hat = wmean(y.*sinwt_tau,wei);
    
    YC = YC_hat - Y*C;
    YS = YS_hat - Y*S;
    
    CC_hat = wmean(coswt_tau.^2,wei);
    SS_hat = wmean(sinwt_tau.^2,wei);
    
    CC = CC_hat - C^2;
    SS = SS_hat - S^2;
    
    gls(ifreq) = 1/YY*(YC^2/CC+YS^2/SS);
    
end

end