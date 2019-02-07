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
function s = mwi_model_2T13T2scc_epgx_moreoption(fa,te,tr,...
                                      Amw,Aiw,Aew,...
                                      t2smw,t2siw,t2sew,...
                                      t1s,t1l,...
                                      fmw,fiw,few,...
                                      totalfield,pini,b1,ka,npulse,RFphi0,T3D_all,...
                                      model)
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
if nargin < 18
%     ka = 2e-3;    %ms^-1
    ka=2;   % s^-1
end
if nargin < 19
    npulse = 200;
end
if nargin < 20
    RFphi0 = 50;
end
if nargin < 22
    model = 'epgx';
end


%% T1 weighting with RF spoiling
% EPG-X
% phi0 = 50;      % initial RF phase
phiCycle = RF_phase_cycle(npulse,RFphi0);
t1x = [t1l, t1s]; 
t2x = [t2siw,t2smw]; % assuming T2* of iw has similar T2 of long T1 compartment
fx = Amw/(Aiw+Aew+Amw); % mwf
fs = (fmw-fiw); % frequency difference between long and short T1 compartments
% fs = (fmw-(fiw*Aiw+few*Aew)/(Aiw+Aew)); % frequency difference between long and short T1 compartments
nfa=length(fa);
SF = zeros(nfa,2);
for ii=1:nfa
    % true flip angle
    alpha = fa(ii)*b1;
   
    % Compute RF spoling phase cycles
    % 2 pools, with exchange
    % start with steady-state signals, transverse magnitisation
    s1 = Signal_GRE_T1wMono((1-fx), alpha, t1x(1), tr);
    s2 = Signal_GRE_T1wMono(fx, alpha, t1x(2), tr);
    
    % Note: minimum effect if abs. is taken here
    % saturation factors
    switch lower(model)
        case 'epgx'
            % EPG-X core
            z1 = s1/sind(alpha);    % long T1 compartment, longitudinal magnitisation (z-state)
            z2 = s2/sind(alpha);    % short T1 compartment, longitudinal magnitisation (z-state)
            tmp = EPGX_GRE_BMsplit_PrecomputedT(T3D_all{ii},phiCycle,tr,t1x,t2x,fx,ka,'delta',fs,'kmax',10,'ss',[z1,z2]);
%             SF(ii,1) = abs(tmp{2}(end));
%             SF(ii,2) = abs(tmp{1}(end)); 
            SF(ii,1) = (tmp{2}(end));
            SF(ii,2) = (tmp{1}(end));
%             SF(ii,1) = conj(tmp{2}(end));
%             SF(ii,2) = conj(tmp{1}(end));

        case 'epg'
            AA = d2r(alpha)*ones(npulse,1);
%             tmp{1} = EPG_GRE(AA,phiCycle,tr,t1x(1),t2x(1),'kmax',10);
%             tmp{2} = EPG_GRE(AA,phiCycle,tr,t1x(2),t2x(2),'kmax',10);
%             z1 = s1/sind(alpha);    % long T1 compartment, longitudinal magnitisation (z-state)
%             z2 = s2/sind(alpha);    % short T1 compartment, longitudinal magnitisation (z-state)
%             tmp2{1} = EPG_GRE_precomputedT(T3D_all{ii},AA,phiCycle,tr,t1x(1),t2x(1),'kmax',10,'ss',z1);
%             tmp2{2} = EPG_GRE_precomputedT(T3D_all{ii},AA,phiCycle,tr,t1x(2),t2x(2),'kmax',10,'ss',z2);
            tmp{1} = EPG_GRE_precomputedT(T3D_all{ii},AA,phiCycle,tr,t1x(1),t2x(1),'kmax',10);
            tmp{2} = EPG_GRE_precomputedT(T3D_all{ii},AA,phiCycle,tr,t1x(2),t2x(2),'kmax',10);
            SF(ii,1) = (tmp{2}(end))*fx;
            SF(ii,2) = (tmp{1}(end))*(1-fx); 
            
        case 'standard'
            SF(ii,1) = s2;
            SF(ii,2) = s1;
    end
     
end

% initiate 2D signal with the saturation factors
[SFmy,teM] = ndgrid(SF(:,1),te);
[SFl,~] = ndgrid(SF(:,2),te);
if length(totalfield) == nfa
    [totalfield,~] = ndgrid(totalfield,te);
end
if length(pini) == nfa
    [pini,~] = ndgrid(pini,te);
end

%% T2* weighting applied here
s = (Amw+Aiw+Aew) * ...
    (SFmy               .*exp(-teM*(1/t2smw-1i*2*pi*fmw)) + ...     % mw
     SFl*(Aiw/(Aiw+Aew)).*exp(-teM*(1/t2siw-1i*2*pi*fiw)) + ...     % iw
     SFl*(Aew/(Aiw+Aew)).*exp(-teM*(1/t2sew-1i*2*pi*few))).* ...    % ew
     exp(1i*2*pi*totalfield.*teM).*exp(-1i*pini);                         % initial phase and total field

end