%% function output = function_name(input)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function [fim, crlb] = mwi_crlb_3cc_mwf_GivenicvfHCFM(te,sigma,S0,mwf,icvf,r2s_mw,r2s_iw,r2s_ew,w_mb,w_ib,w_eb,phi)

fim = zeros(7);
for kt = 1:length(te)
    fim = fim + FIMx(te(kt));
end

crlb = diag(pinv(fim));

    function FIMij = FIMx(t)
        dSdS0   = exp(1i*phi)*(mwf*exp(-r2s_mw*t + 1i*w_mb*t) + ...
                               (1-mwf)*icvf*exp(-r2s_iw*t + 1i*w_ib*t) + ...
                               (1-mwf)*(1-icvf)*exp(-r2s_ew*t + 1i*w_eb*t));
        dSdmwf  = exp(1i*phi)*S0*(exp(-r2s_mw*t + 1i*w_mb*t) - ...
                                      icvf*exp(-r2s_iw*t + 1i*w_ib*t) - ...
                                      (1-icvf)*exp(-r2s_ew*t + 1i*w_eb*t));
%         dSdicvf = exp(1i*phi)*S0*((1-mwf)*exp(-r2s_iw*t + 1i*w_ib*t) - ...
%                                               (1-mwf)*exp(-r2s_ew*t + 1i*w_eb*t));

        dSdr2s_mw = -t * exp(1i*phi)*S0*mwf * exp(-r2s_mw*t + 1i*w_mb*t);
        dSdr2s_iw = -t * exp(1i*phi)*S0*(1-mwf)*icvf * exp(-r2s_iw*t + 1i*w_ib*t);
        dSdr2s_ew = -t * exp(1i*phi)*S0*(1-mwf)*(1-icvf) * exp(-r2s_ew*t + 1i*w_eb*t);

%         dSdw_mb = 1i*t * exp(1i*phi)*S0*mwf * exp(-r2s_mw*t + 1i*w_mb*t);
%         dSdw_ib = 1i*t * exp(1i*phi)*S0*(1-mwf)*icvf * exp(-r2s_iw*t + 1i*w_ib*t);
        dSdw_eb = 1i*t * exp(1i*phi)*S0*(1-mwf)*(1-icvf) * exp(-r2s_ew*t + 1i*w_eb*t);

        dSdphi = 1i*exp(1i*phi)*S0*(mwf*exp(-r2s_mw*t + 1i*w_mb*t) + ...
                                    (1-mwf)*icvf*exp(-r2s_iw*t + 1i*w_ib*t) + ...
                                    (1-mwf)*(1-icvf)*exp(-r2s_ew*t + 1i*w_eb*t));

%         dSdtheta = [dSdS0;dSdmwf;dSdicvf;dSdr2s_mw;dSdr2s_iw;dSdr2s_ew;dSdw_mb;dSdw_ib;dSdw_eb;dSdphi];
    dSdtheta = [dSdS0;dSdmwf;dSdr2s_mw;dSdr2s_iw;dSdr2s_ew;dSdw_eb;dSdphi];

        FIMij = (1/(sigma.^2))*real(dSdtheta * dSdtheta');
        
    end
end

