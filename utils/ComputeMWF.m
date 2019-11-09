%% function mwf = ComputeMWF(fitRes)
%
% Input
% --------------
% fitRes        : struct array from mwi process
%
% Output
% --------------
% mwf           : myelin water fraction
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 3 August 2018
% Date last modified: 25 October 2019
%
%
function [mwf,error_map] = ComputeMWF(fitRes)

if isfield(fitRes,'S0_MW')
    if isfield(fitRes,'S0_IEW')
        S0 = (fitRes.S0_MW + fitRes.S0_IEW);
    else
        S0 = (fitRes.S0_MW + fitRes.S0_IW + fitRes.S0_EW);
    end
    mwf = fitRes.S0_MW ./ S0;
else
    if size(fitRes.estimates,4) < 8
        S0 = sum(fitRes.estimates(:,:,:,1:2),4);
    else
        S0 = sum(fitRes.estimates(:,:,:,1:3),4);
    end
    mwf = (fitRes.estimates(:,:,:,1)) ./ S0;
end

nanMap = isnan(mwf) ;
infMap = isinf(mwf) ;
nullSignal = S0 == 0;
error_map = or(nanMap,infMap) .* nullSignal;

mwf(isnan(mwf)) = 0;
mwf(isinf(mwf)) = 0;

end