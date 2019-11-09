%% [fim, crlb] = mwi_crlb_3cc_mwf_GivenicvfHCFM(te,sigma,S0,mwf,v_ic,r2s_mw,r2s_iw,r2s_ew,w_mb,w_ib,w_eb,phi)
%
% Input
% --------------
% te            : echo times, in second
% sigma         : noise SD, scalar or vector
% S0            : signal intensity of all compartments at TE=0
% mwf           : myelin water fraction
% v_ic          : volume fraction of intra-axonal water (=Aiw/(Aiw+Aew))
% r2s_mw        : R2* of myelin water, in s^-1
% r2s_iw        : R2* of intra-axonal water, in s^-1
% r2s_ew        : R2* of extracellular water, in s^-1
% w_mb          : angular frequency shift of myelin water with background,
%                 in rad (=2*pi*(fmw+fbkg))
% w_ib          : angular frequency shift of intra-axonal water with background,
%                 in rad (=2*pi*(fiw+fbkg))
% w_eb          : angular frequency shift of extracellular water with background,
%                 in rad (=2*pi*(few+fbkg))
% phi           : initial phase offset
%
% Output
% --------------
% fim           : Fisher information matrix
% crlb          : Cramer-rao lower bound of each variables
%
% Description: Compute the Cramer-rao lower bound of 3-pool myelin water
% imaging (MWI) model with complex-valued fitting
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 17 October 2019
% Date last modified:
%
%
function [fim, crlb] = mwi_crlb_3cc_mwf_GivenicvfHCFM(te,sigma,S0,mwf,v_ic,r2s_mw,r2s_iw,r2s_ew,w_mw,w_iw,w_ew,w_bkg,phi)

% if sigma is a scalar, that means the noise is the same across echoes
if isscalar(sigma)
    sigma = ones(size(te)) * sigma;
end

% combine background field and compartmental shift
w_mb = w_mw + w_bkg;
w_ib = w_iw + w_bkg;
w_eb = w_ew + w_bkg;

% compute the Fisher information matrix
fim = zeros(7);
for kt = 1:length(te)
    fim = fim + FIMx(te(kt), sigma(kt));
end

% compute the inverse of Fisher information matrix
crlb = diag(pinv(fim));

    % nested function to compute Fisher information matrix
    function FIMij = FIMx(t, sigma)
        % partial derivative w.r.t. each variable
        dSdS0   = exp(1i*phi) * ( mwf     *            exp(-r2s_mw*t + 1i*w_mb*t) + ...
                                  (1-mwf) * v_ic     * exp(-r2s_iw*t + 1i*w_ib*t) + ...
                                  (1-mwf) * (1-v_ic) * exp(-r2s_ew*t + 1i*w_eb*t));
        dSdmwf  = exp(1i*phi)*S0 * (           exp(-r2s_mw*t + 1i*w_mb*t) - ...
                                    v_ic     * exp(-r2s_iw*t + 1i*w_ib*t) - ...
                                    (1-v_ic) * exp(-r2s_ew*t + 1i*w_eb*t));
%         dSdvic = exp(1i*phi)*S0*(1-mwf) * (exp(-r2s_iw*t + 1i*w_ib*t) - ...
%                                            exp(-r2s_ew*t + 1i*w_eb*t));

        dSdr2s_mw = -t * exp(1i*phi)*S0 * mwf                * exp(-r2s_mw*t + 1i*w_mb*t);
        dSdr2s_iw = -t * exp(1i*phi)*S0 * (1-mwf) * v_ic     * exp(-r2s_iw*t + 1i*w_ib*t);
        dSdr2s_ew = -t * exp(1i*phi)*S0 * (1-mwf) * (1-v_ic) * exp(-r2s_ew*t + 1i*w_eb*t);

%         dSdw_mb = 1i*t * exp(1i*phi)*S0 * mwf                * exp(-r2s_mw*t + 1i*w_mb*t);
%         dSdw_ib = 1i*t * exp(1i*phi)*S0 * (1-mwf) * v_ic     * exp(-r2s_iw*t + 1i*w_ib*t);
%         dSdw_eb = 1i*t * exp(1i*phi)*S0 * (1-mwf) * (1-v_ic) * exp(-r2s_ew*t + 1i*w_eb*t);
        
        dSdw_bkg = 1i*t*exp(1i*phi)*S0 * ( mwf               * exp(-r2s_mw*t + 1i*w_mb*t) + ...
                                          (1-mwf) * v_ic     * exp(-r2s_iw*t + 1i*w_ib*t) + ...
                                          (1-mwf) * (1-v_ic) * exp(-r2s_ew*t + 1i*w_eb*t));

        dSdphi = 1i*exp(1i*phi)*S0 * ( mwf                * exp(-r2s_mw*t + 1i*w_mb*t) + ...
                                       (1-mwf) * v_ic     * exp(-r2s_iw*t + 1i*w_ib*t) + ...
                                       (1-mwf) * (1-v_ic) * exp(-r2s_ew*t + 1i*w_eb*t));

%         dSdtheta = [dSdS0;dSdmwf;dSdvic;dSdr2s_mw;dSdr2s_iw;dSdr2s_ew;dSdw_mb;dSdw_ib;dSdw_eb;dSdphi];
%         dSdtheta = [dSdS0;dSdmwf;dSdr2s_mw;dSdr2s_iw;dSdr2s_ew;dSdw_eb;dSdphi];
        dSdtheta = [dSdS0;dSdmwf;dSdr2s_mw;dSdr2s_iw;dSdr2s_ew;dSdw_bkg;dSdphi];

        FIMij = (1/(sigma.^2))*real(dSdtheta * dSdtheta');
        
    end
end

