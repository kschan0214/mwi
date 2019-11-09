%% s = mwi_model_3cc_dimwi(te,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI)
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
% Results (perp./para.) using 16 echoes:
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
function s = mwi_model_3cc_dimwi(te,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI)

%% validate input variables
if nargin < 13
	DIMWI.isFreqMW  = false;
    DIMWI.isFreqIW  = false;
    DIMWI.isR2sEW   = false;
end
if nargin < 12
    pini=0;
end

if DIMWI.isFreqMW || DIMWI.isFreqIW || DIMWI.isR2sEW
    fixParam = DIMWI.fixParam;
end

% combine compartmental frequency with background
freq_mwbg = freq_mw + fbg;
freq_iwbg = freq_iw + fbg;
freq_ewbg = freq_ew + fbg;

%% check DIMWI model
%%%%%%%%%% frequency shifts estimated using HCFM %%%%%%%%%%
if DIMWI.isFreqMW || DIMWI.isFreqIW
    % derive g-ratio 
    g = sqrt(Aiw/(Aiw+Amw/fixParam.rho_mw));
    % compute frequency shift given theta
    [freq_mw, freq_iw] = hcfm_freq(fixParam.x_i,fixParam.x_a,g,DIMWI.theta,fixParam.E,fixParam.b0);
    if DIMWI.isFreqMW 
        freq_mwbg = freq_mw + fbg;
    end
    if DIMWI.isFreqMW 
        freq_iwbg = freq_iw + fbg;
    end
end

%%%%%%%%%% extra decay on extracellular water estimated by HCFM %%%%%%%%%%
if DIMWI.isR2sEW
    
    % assume extracellular water has the same T2* as intra-axonal water
    t2sew = t2siw;
    
    % derive fvf from signal intensity
    v_ic = Aiw ./ (Aiw + Aew);
    fvf = v_ic ./ (g^2 - v_ic*g^2 + v_ic); 
    
    % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
    d_e = hcfm_dephasing_decay_DE(te,fvf,g,fixParam.x_i,fixParam.x_a,DIMWI.theta,fixParam.b0);
    
else
    d_e = 0;
end

%%

s = (Amw*exp(te*(-1/t2smw+1i*2*pi*freq_mwbg)) + ...
     Aiw*exp(te*(-1/t2siw+1i*2*pi*freq_iwbg)) + ...
     Aew*exp(te*(-1/t2sew+1i*2*pi*freq_ewbg)).*exp(-d_e))*exp(1i*pini);

end