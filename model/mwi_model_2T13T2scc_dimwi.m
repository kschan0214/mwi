%% s = mwi_model_2T13T2scc_dimwi(te,tr,fa,b1,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1mw,t1iew,kiewmw,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI,EPGX)
%
% Input (initial guess,[lb,ub])
% --------------
% te            : echo times, in second
% Amw           : Myelin water signal
% Aiw           : Axonal water signal
% Aew           : extracellular water signal
% t2smw         : myelin water T2*, in second
% t2siw         : axonal water T2*, in second
% t2sew         : extracellular water T2*, in second
% freq_mw    	: myelin water frequency, in Hz
% freq_iw   	: axonal water frequency, in Hz
% freq_ew       : extracellular frequency, in Hz
% fbg           : background field, in Hz
% pini          : initial phase introduced by B1+ phase offset, in rad
% DIMWI         : structure for diffusion informed MWI
%
% Output
% --------------
% s             : complex-valued 3-pool signal
%
% Description: Complex-valued model fitted to complex-valued data used in
% Nam et al. NeuroImages 2015 116:214-221
% Protocol (3T):
%   - voxel size    = 2mm isotropic
%   - fa            = 30
%   - TR            = 120ms
%   - TE            = 2.1:1.93:61.93ms (32 echoes)
%   - BW            = 1502 Hz/pixel
% Key processing:
%   (1) averaged two adjacent slices, weighted ls fitting with weight of
%   magn (voxel-by-voxel)
%   (2) averaged magn., angle of ROI-averaged complex data (ROI)
% Resu created: 4 January 2018
% Date last modified: 16 August 2018
% Date last modified: 29 October 2019
%
%lts (perp./para.) using 16 echoes:
%   - MWF           = 14.4/8.5 %
%   - f_(my-ex)     = 11.1/2.5 Hz
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 4 January 2018
% Date last modified: 16 August 2018
% Date last modified: 29 October 2019
%
%
function s = mwi_model_2T13T2scc_dimwi(te,tr,fa,b1,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1mw,t1iew,kiewmw,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI,EPGX)

%% validate input variables
if nargin < 20
    EPGX.isExchange = false;
    EPGX.isEPG = false;
end
    
if nargin < 19 || isempty(DIMWI)
	DIMWI.isFreqMW  = false;
    DIMWI.isFreqIW  = false;
    DIMWI.isR2sEW   = false;
end
if nargin < 18
    pini=0;
end

if DIMWI.isFreqMW || DIMWI.isFreqIW || DIMWI.isR2sEW
fixParam = DIMWI.fixParam;
end

if numel(pini) == 1
    pini = pini * ones(size(fa));
end
if numel(fbg) == 1
    fbg = fbg * ones(size(fa));
end

%% check DIMWI model
%%%%%%%%%% frequency shifts estimated using HCFM %%%%%%%%%%
if DIMWI.isFreqMW || DIMWI.isFreqIW
    % derive g-ratio 
    g = sqrt(abs(Aiw)/(abs(Aiw)+abs(Amw)/fixParam.rho_mw));
    % compute frequency shift given theta
    [freqMW, freqIW] = hcfm_freq(fixParam.x_i,fixParam.x_a,g,DIMWI.theta,fixParam.E,fixParam.b0);
    if DIMWI.isFreqMW 
        freq_mw = freqMW;
    end
    if DIMWI.isFreqMW 
        freq_iw = freqIW;
    end
end

%%%%%%%%%% extra decay on extracellular water estimated by HCFM %%%%%%%%%%
if DIMWI.isR2sEW
    
    % assume extracellular water has the same T2* as intra-axonal water
    t2sew = t2siw;
    
    % derive fvf from signal intensity
    v_ic = abs(Aiw) ./ (abs(Aiw) + abs(Aew));
    fvf = v_ic ./ (g^2 - v_ic*g^2 + v_ic); 
    
    % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
    d_e = hcfm_dephasing_decay_DE(te,fvf,g,fixParam.x_i,fixParam.x_a,DIMWI.theta,fixParam.b0);
    
else
    d_e = zeros(size(te));
end

%% T1
if EPGX.isExchange
    if EPGX.isEPG
        % EPG-X
        phiCycle = RF_phase_cycle(EPGX.npulse,EPGX.rfphase);
        t1x = [t1iew, t1mw]; 
        t2x = [t2siw,t2smw]; % assuming T2* of iw has similar T2 of long T1 compartment
%         t2x = [t2siw,284e-6]; % assuming T2* of iw has similar T2 of long T1 compartment
        fx = (Amw/EPGX.rho_mw)/(Aiw+Aew+Amw/EPGX.rho_mw); % mwf
        fs = (freq_mw-freq_iw); % frequency difference between long and short T1 compartments
        
        SF = zeros(length(fa),2);
        for ii=1:length(fa)
            % true flip angle
            alpha = fa(ii)*b1;

            % 2 pools, with exchange
            % start with steady-state signals, longitudinal magnetisation (Mz)
            z1 = Signal_GRE_T1wMono(1-fx, alpha, t1iew, tr)/sind(alpha);
            z2 = Signal_GRE_T1wMono(fx,	alpha, t1mw, tr)/sind(alpha);
            % EPG-X core
            tmp = EPGX_GRE_BMsplit_PrecomputedT(EPGX.T3D_all{ii},phiCycle,tr,t1x,t2x,fx,kiewmw,'delta',fs,'kmax',10,'ss',[z1,z2]);

            % saturation factors  
            SF(ii,1) = tmp{2}(end);
            SF(ii,2) = tmp{1}(end); 
        end
        ss(1,:) = SF(:,1) * (Aiw+Aew+Amw) * EPGX.rho_mw;
        ss(2,:) = SF(:,2) * (Aiw+Aew+Amw) * (Aiw/(Aiw+Aew));
        ss(3,:) = SF(:,2) * (Aiw+Aew+Amw) * (Aew/(Aiw+Aew));
    else
        % Exchange only
        [~, ss] = mwi_model_ssSPGR_2T1(fa*b1,tr,Amw/0.43,Aiw+Aew,t1mw,t1iew,[],kiewmw);
        ss(1,:) = ss(1,:)* EPGX.rho_mw;
        ss(3,:) = Aew/(Aiw+Aew) * ss(2,:);
        ss(2,:) = Aiw/(Aiw+Aew) * ss(2,:);
    end
else
    % separated steady-state
    ss = mwi_model_3T1_ssSPGR(fa,tr,Amw,Aiw,Aew,t1mw,t1iew,t1iew,b1);
end

[te2D,Amw2D]    = ndgrid(te,ss(1,:));
[~,Aiw2D]       = ndgrid(te,ss(2,:));
[~,Aew2D]       = ndgrid(te,ss(3,:));
[~,fbg2D]       = ndgrid(te,fbg);
[~,pini2D]      = ndgrid(te,pini);
[d_e2D,~]       = ndgrid(d_e,fa);


%%

s = (Amw2D.*exp(te2D*(-1/t2smw+1i*2*pi*freq_mw)) + ...
     Aiw2D.*exp(te2D*(-1/t2siw+1i*2*pi*freq_iw)) + ...
     Aew2D.*exp(te2D*(-1/t2sew+1i*2*pi*freq_ew)).*exp(-d_e2D)).*exp(1i*pini2D).*exp(1i*2*pi*fbg2D.*te2D);

end