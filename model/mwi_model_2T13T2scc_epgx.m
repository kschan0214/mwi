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
function s = mwi_model_2T13T2scc_epgx(fa,te,tr,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1s,t1l,fmw,fiw,few,pini,b1,ka,npulse,T3D_all)
% set default values
if nargin < 15
    pini=0;
end
if nargin < 16
    b1=1;
end
if nargin < 17
%     ka = 2e-3;    %ms^-1
    ka=2;   % s^-1
end
if nargin < 18
    npulse = 200;
end

%% T1 weighting with RF spoiling
% EPG-X
phi0 = 50;      % initial RF phase
phiCycle = RF_phase_cycle(npulse,phi0);
t1x = [t1l, t1s]; 
t2x = [t2siw,t2smw]; % assuming T2* of iw has similar T2 of long T1 compartment
fx = Amw/(Aiw+Aew+Amw); % mwf
fs = (fmw-fiw); % frequency difference between long and short T1 compartments

nfa=length(fa);
SF = zeros(nfa,2);
for ii=1:nfa
    % true flip angle
    alpha = fa(ii)*b1;

   
    % Compute RF spoling phase cycles
    % 2 pools, with exchange
    % start with steady-state signals
    z1 = Signal_GRE_T1wMono((1-fx), alpha, t1x(1), tr)/sind(alpha);
    z2 = Signal_GRE_T1wMono(fx, alpha, t1x(2), tr)/sind(alpha);
    % EPG-X core
    tmp = EPGX_GRE_BMsplit_PrecomputedT(T3D_all{ii},phiCycle,tr,t1x,t2x,fx,ka,'delta',fs,'kmax',10,'ss',[z1,z2]);
    
    % saturation factors
    % Note: no effect if abs is taken here
    SF(ii,1) = abs(tmp{2}(end));
    SF(ii,2) = abs(tmp{1}(end));    
end

% initiate 2D signal with the saturation factors
[SFmy,teM] = ndgrid(SF(:,1),te);
[SFl,~] = ndgrid(SF(:,2),te);

%% T2* weighting applied here
s = (Amw+Aiw+Aew) * ...
    (SFmy               .*exp(-teM*(1/t2smw-1i*2*pi*fmw)) + ... % mw
     SFl*(Aiw/(Aiw+Aew)).*exp(-teM*(1/t2siw-1i*2*pi*fiw)) + ... % iw
     SFl*(Aew/(Aiw+Aew)).*exp(-teM*(1/t2sew-1i*2*pi*few)))* ... % ew
     exp(-1i*pini);                                             % initial phase

end