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
function s = mwi_model_2cm_jointT1T2s_epg(fa,te,tr,As,Al,t2ss,t2sl,t1s,t1l,fs,fl,b1)

if nargin < 15
    b1=1;
end

% EPG
phiCycle = 50;
nfa=length(fa);
npulse = floor(5*t1l/tr);
phiCycle = RF_phase_cycle(npulse,phiCycle);
t2s = 15e-3;
t2l = 100e-3;

SF = zeros(nfa,2);
for ii=1:nfa
    alpha = fa(ii)*b1;
    AA = d2r(alpha)*ones(npulse,1);
   
    % Compute RF spoling phase cycles
    % single pool, no exchange
    tmp1 = EPG_GRE(AA,phiCycle,tr,t1s,t2s);
    tmp2 = EPG_GRE(AA,phiCycle,tr,t1l,t2l);
    
    % Note: no effect if abs is taken
    SF(ii,1) = abs(tmp1(end));
    SF(ii,2) = abs(tmp2(end));    
end

% T1w signal with RF spoiling effect
[SFs,teM] = ndgrid(SF(:,1),te);
[SFl,~] = ndgrid(SF(:,2),te);

% magnitude combined signal
s = abs(SFs*As.*exp(-teM*(1/t2ss+1i*2*pi*fs)) + SFl*Al.*exp(-te*(1/t2sl+1i*2*pi*fl)));


end