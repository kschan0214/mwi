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
% Date last modified:
%
%
function [mwf,wf] = ComputeMWF(fitRes)

mwf = (fitRes.estimates(:,:,:,1)) ./ sum(fitRes.estimates(:,:,:,1:3),4);

mwf(isnan(mwf)) = 0;
mwf(isinf(mwf)) = 0;

iwf = (fitRes.estimates(:,:,:,2)) ./ sum(fitRes.estimates(:,:,:,1:3),4);
ewf = (fitRes.estimates(:,:,:,3)) ./ sum(fitRes.estimates(:,:,:,1:3),4);
wf.mwf = mwf;
wf.iwf = iwf;
wf.ewf = ewf;

end