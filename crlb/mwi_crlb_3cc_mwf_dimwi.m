%% [fim, crlb] = mwi_crlb_3cc_mwf_dimwi(te,sigma,S0,mwf,v_ic,theta,r2s_mw,r2s_iw,freq_ew,freq_bkg,phi,b0,fixParam)
%
% Input
% --------------
% te            : echo times, in second
% sigma         : noise SD, scalar or vector
% S0            : sum of signal intensity of all compartments at TE=0
% mwf           : myelin water fraction
% v_ic          : volume fraction of intra-axonal water (=Aiw/(Aiw+Aew))
% theta         : angle between B0 and fibre orientation, in rad
% r2s_mw        : R2* of myelin water, in s^-1
% r2s_iw        : R2* of intra-axonal water, in s^-1
% freq_ew       : frequency shift of extracellular water, in Hz
% freq_bkg     	: frequency shift of background, in Hz
% phi           : initial phase offset, in rad
% b0            : field strength, in Tesla
% fixParam      : structure contain the following constants
%   rho_mw      : relative myelin water density
%   x_i         : myelin shealth isotropic susceptibility, in ppm
%   x_a         : myelin shealth anisotropic susceptibility, in ppm
%   E           : exchange term, in ppm
%
% Output
% --------------
% fim           : Fisher information matrix
% crlb          : Cramer-rao lower bound of each variables
%
% Description: Compute the Cramer-rao lower bound of 3-pool diffusion informed myelin water
% imaging (DIMWI) model with complex-valued fitting
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 29 October 2019
% Date last modified:
%
%
function [fim, crlb] = mwi_crlb_3cc_mwf_dimwi(te,sigma,S0,mwf,v_ic,theta,r2s_mw,r2s_iw,freq_ew,freq_bkg,phi,b0,fixParam)

%%%%%%%%%% check input %%%%%%%%%%
% fixed parameters
try rho_mw  = fixParam.rho_mw;  catch; rho_mw   = 0.43; end
try x_i     = fixParam.x_i;     catch; x_i      = -0.1; end
try x_a     = fixParam.x_a;     catch; x_a      = -0.1; end
try E       = fixParam.E;       catch; E        = 0.02; end

% if sigma is a scalar, that means the noise is the same across echoes
if isscalar(sigma)
    sigma = ones(size(te)) * sigma;
end

%%%%%%%%%% convert parameters in correct convention %%%%%%%%%%
% derive g-ratio and fvf based on input
g   = sqrt((1-mwf)*v_ic/((1-mwf)*v_ic+mwf/rho_mw));
fvf = v_ic / (g^2 - v_ic*g^2 + v_ic); 

% signal dephase in extracellular water due to myelin sheath, Eq.[A7]
D_E = hcfm_dephasing_decay_DE(te,fvf,g,x_i,x_a,theta,b0);

% compute frequency shift given theta
[freq_mw, freq_iw] = hcfm_freq(x_i,x_a,g,theta,E,b0);
% combine background field and compartmental shift
w_mw    = freq_mw*2*pi;
w_iw    = freq_iw*2*pi;
w_ew    = freq_ew*2*pi;
w_bkg   = freq_bkg*2*pi;

%%%%%%%%%% CRLB main %%%%%%%%%%
% compute the Fisher information matrix
fim = zeros(6);
for kt = 1:length(te)
    fim = fim + FIMx(te(kt), sigma(kt), D_E(kt));
end

% compute the inverse of Fisher information matrix
crlb = diag(pinv(fim));

    % nested function to compute Fisher information matrix
    function FIMij = FIMx(t, sigma, D_Et)
        %%%%%%%%%% partial derivative w.r.t. each variable %%%%%%%%%%
        % Signal compartment 
        dSdS0   = exp(1i*phi)*exp(1i*w_bkg*t) * ...                         % external
            ( mwf     *            exp(-r2s_mw*t + 1i*w_mw*t) + ...         % myelin water
              (1-mwf) * v_ic     * exp(-r2s_iw*t + 1i*w_iw*t) + ...         % intra-axonal water
              (1-mwf) * (1-v_ic) * exp(-r2s_iw*t + 1i*w_ew*t)*exp(-D_Et));  % extracellular water
        
        
        dSdmwf  = S0*exp(1i*phi)*exp(1i*w_bkg*t) * ...
                                    (           exp(-r2s_mw*t + 1i*w_mw*t) - ...
                                     v_ic     * exp(-r2s_iw*t + 1i*w_iw*t) - ...
                                     (1-v_ic) * exp(-r2s_iw*t + 1i*w_ew*t)*exp(-D_Et));

%         dSdv_ic = S0*exp(1i*phi)*exp(1i*w_bkg*t) * (1-mwf)*exp(-r2s_iw*t)*...
%                                     (exp(1i*w_iw*t) - exp(1i*w_ew*t)*exp(-D_Et));
        
        % R2*
        dSdr2s_mw = S0*exp(1i*phi)*exp(1i*w_bkg*t) * ...
                        (-t * mwf * exp(-r2s_mw*t + 1i*w_mw*t));
                    
        dSdr2s_iw = S0*exp(1i*phi)*exp(1i*w_bkg*t) * (-t)*(1-mwf)*exp(-r2s_iw*t) * ...
                        (v_ic * exp(1i*w_iw*t) + (1-v_ic)*exp(1i*w_ew*t)*exp(-D_Et));
        
%         % frequency shift
%         dSdw_mw = -1i*dSdr2s_mw;
%         dSdw_iw = S0*exp(1i*phi)*exp(1i*w_bkg*t) * (1i*t*(1-mwf)*v_ic*exp(-r2s_iw*t + 1i*w_iw*t));
%         dSdw_ew = S0*exp(1i*phi)*exp(1i*w_bkg*t) * (1i*t*(1-mwf)*(1-v_ic)*exp(-r2s_iw*t + 1i*w_ew*t)*exp(-D_Et));
        
        % external 
        dSdphi      = 1i*S0*dSdS0;
        dSdw_bkg    = t*dSdphi;

        % unknown variables vector
        dSdtheta = [dSdS0;dSdmwf;dSdr2s_mw;dSdr2s_iw;dSdw_bkg;dSdphi];

        FIMij = (1/(sigma.^2))*real(dSdtheta * dSdtheta');
        
    end
end

