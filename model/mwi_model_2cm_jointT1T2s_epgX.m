%% s = mwi_model_2cm_jointT1T2s_epg(fa,te,tr,As,Al,t2ss,t2sl,t1s,t1l,fs,fl,b1)
%
% Input
% --------------
% fa            : flip angles
% te            : echo times
% tr            : repetition time
% Amy           : Myelin water signal
% Aax           : Axonal water signal
% Aex           : extracellular water signal
% t2smy         : myelin water T2*
% t2sax         : axonal water T2*
% t2sex         : extracellular water T2*
% t1my          : myelin water T1
% t1ax          : axonal water T1
% t1ex          : extracellular water T1
% fmy           : myelin water frequency - extracellular water
% fax           : axonal water frequency - extracellular water
% b1            : factor of B1+ inhomogeneity: true flip angle=1
%
% Output
% --------------
% s             : 2D magnitude 3-pool signal, 1st dim=fa;2nd dim=te
%
% Description: extended magnitude signal model account for the effects
% of T1w
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 18 January 2018
% Date last modified:
%
%
% function s = mwi_model_2cm_jointT1T2s_epgX(fa,te,tr,As,Al,t2ss,t2sl,t1s,t1l,fs,fl,b1,ka)
function s = mwi_model_2cm_jointT1T2s_epgX(fa,te,tr,A0,mwf,t2ss,t2sl,t1s,t1l,fs,fl,b1,ka)

if nargin < 12
    b1=1;
end
if nargin < 13
%     kx = 2e-3;    %ms^-1
    ka=2;   % s^-1
end

%% T1 weighting with RF spoiling
% EPG
phi0 = 50;      % initial RF phase
nfa=length(fa);
% npulse = floor(5*t1l/tr);
npulse = 150;  % signal usually reaches SS quickly(~80), preset the no.of pulses to avoid long computational time
phiCycle = RF_phase_cycle(npulse,phi0);
t1x = [t1l t1s];
t2x = [100 25]*1e-3;    
% fx = As/(Al+As);
fx = mwf;

SF = zeros(nfa,2);
for ii=1:nfa
    alpha = fa(ii)*b1;
    AA = d2r(alpha)*ones(npulse,1);
   
    % Compute RF spoling phase cycles
    % 2 pools, exchange
    tmp = EPGX_GRE_BMsplit(AA,phiCycle,tr,t1x,t2x,fx,ka,'delta',fs);
    
    % Note: no effect if abs is taken here
    SF(ii,1) = abs(tmp{2}(end));
    SF(ii,2) = abs(tmp{1}(end));    
    
end

% T1w signal with RF spoiling effect
[SFs,teM] = ndgrid(SF(:,1),te);
[SFl,~] = ndgrid(SF(:,2),te);

% magnitude combined signal
% s = abs(SFs*As.*exp(-teM*(1/t2ss+1i*2*pi*fs)) + SFl*Al.*exp(-te*(1/t2sl+1i*2*pi*fl)));

% s = (As+Al)*abs(SFs.*exp(-teM*(1/t2ss+1i*2*pi*fs)) + SFl.*exp(-teM*(1/t2sl+1i*2*pi*fl)));
s = A0*abs(SFs.*exp(-teM*(1/t2ss+1i*2*pi*fs)) + SFl.*exp(-teM*(1/t2sl+1i*2*pi*fl)));


end