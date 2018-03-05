%% B = mwi_GetPhaseAndDecay(t,freqs,r2s)
%
% Input
% --------------
% fa            : flip angle (1xm)
% tr            : repetition time
% t1            : species T1 (1xn)       
%
% Output
% --------------
% B             : mxn matrix signal T1 decay across flip angle
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 2 March 2018
% Date last modified:
%
%
function B = mwi_GetSpeciesT1Decay(fa,tr,t1)

% separated species phase and decay
E1 = exp(-tr./t1);

B = sind(fa(:))*(1-E1(:).')./(1-cosd(fa(:))*E1(:).');

end