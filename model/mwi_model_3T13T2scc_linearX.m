%%  mwi_model_2T13T2scc_epgx(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1s,t1l,fmy,fax,b1,ka,npulse,T3D_all)
%
% Input
% --------------
% fa            : flip angles
% te            : echo times
% tr            : repetition time
% Amw           : Myelin water signal
% Aiw           : Intra-axonal water signal
% Aew           : Extracellular water signal
% t2smw         : myelin water T2*
% t2siw         : intra-axonal water T2*
% t2sew         : extracellular water T2*
% t1s           : short T1 (assumed to be myelin water T1) 
% t1l           : long T1 (assumed to be free water (iw & ew) T1)
% fmw           : myelin water + bkg frequency 
% fiw           : intra-axonal water + bkg frequency
% few           : extracellular water + bkg frequency
% b1            : factor of B1+ inhomogeneity: true flip angle=1
% ka            : exchange rate from long T1 to short T1
% npulse        : no. of pulses for EPG-X to reash steady-state
% T3D_all       : transition matrix for EPG-X
%
% Output
% --------------
% s             : 2D complex 3-pool signal, 1st dim=fa;2nd dim=te
%
% Description: extended magnitude signal model account for the effects
% of T1w
% T1,T2* and ka must have the same unit (default assume in second)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 march 2018
% Date last modified: 
%
%
function s = mwi_model_3T13T2scc_linearX(fa,te,tr,...
                                      Amw,Aiw,Aew,Ass,...
                                      t2smw,t2siw,t2sew,...
                                      t1ss,t1w,...
                                      kowmw,kmwss,...
                                      freqmw,freqiw,...
                                      totalfield,pini,b1)
% set default values
if nargin < 15
    totalfield = 0;
end
if nargin < 16
    pini=0;
end
if nargin < 17
    b1=1;
end


fmw = Amw/(Amw+Aiw+Aew+Ass);
fss = Ass/(Amw+Aiw+Aew+Ass);
%% steady state
nfa = length(fa);

alpha = d2r(fa * b1);

M0ow = (1-fmw-fss);
M0mw = fmw;
M0ss = fss;
kmwow = kowmw * M0ow/M0mw;
kssmw = kmwss * M0mw/M0ss;
R1w = 1/t1w;
R1ss = 1/t1ss;

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [[-R1w-kowmw kmwow 0];[kowmw -R1w-kmwow-kmwss kssmw];[0 kmwss -R1ss-kssmw]];
Xi_L = expm(tr*Lambda_L);

% Use same notation as paper
I = eye(3);
C_L = [R1w*M0ow;R1w*M0mw;R1ss*M0ss];
    
SF = zeros(nfa,3);
for kfa = 1:length(fa)
    
    % Use same notation as paper
    T = cos(alpha(kfa))*I;
    
    SF(kfa,:) = sin(alpha(kfa))* inv(I-Xi_L*T) * (Xi_L - I) * inv(Lambda_L) * C_L;

end

[SFow,teM] = ndgrid(SF(:,1),te);
[SFmw,~] = ndgrid(SF(:,2),te);
if length(totalfield) == nfa
    [totalfield,~] = ndgrid(totalfield,te);
end
if length(pini) == nfa
    [pini,~] = ndgrid(pini,te);
end

%% T2* weighting applied here
s = (Amw+Aiw+Aew+Ass) * ...
    (SFmw                .*exp(teM*(-1/t2smw+1i*2*pi*freqmw)) + ...     % mw
     SFow*(Aiw/(Aiw+Aew)).*exp(teM*(-1/t2siw+1i*2*pi*freqiw)) + ...     % iw
     SFow*(Aew/(Aiw+Aew)).*exp(teM*(-1/t2sew))).* ...    % ew
     exp(1i*2*pi*totalfield.*teM).*exp(-1i*pini);                         % initial phase and total field

end