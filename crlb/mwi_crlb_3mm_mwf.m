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
function [fim, crlb] = mwi_crlb_3mm_mwf(te,sigma,S0,mwf,icvf,r2s_mw,r2s_iw,r2s_ew)

fim = zeros(6);
for kt = 1:length(te)
    fim = fim + FIMx(te(kt));
end

crlb = diag(pinv(fim));

    function FIMij = FIMx(t)
        dSdS0   = (mwf*exp(-r2s_mw*t) + ...
                   (1-mwf)*icvf*exp(-r2s_iw*t) + ...
                   (1-mwf)*(1-icvf)*exp(-r2s_ew*t));
        dSdmwf  = S0*(exp(-r2s_mw*t) - ...
                      icvf*exp(-r2s_iw*t) - ...
                      (1-icvf)*exp(-r2s_ew*t));
        dSdicvf = S0*((1-mwf)*exp(-r2s_iw*t) - ...
                      (1-mwf)*exp(-r2s_ew*t));

        dSdr2s_mw = -t * S0*mwf * exp(-r2s_mw*t);
        dSdr2s_iw = -t * S0*(1-mwf)*icvf * exp(-r2s_iw*t);
        dSdr2s_ew = -t * S0*(1-mwf)*(1-icvf) * exp(-r2s_ew*t);

        dSdtheta = [dSdS0;dSdmwf;dSdicvf;dSdr2s_mw;dSdr2s_iw;dSdr2s_ew];

        FIMij = (1/(sigma.^2))*real(dSdtheta * dSdtheta');
        
    end
end

