%% B = mwi_GetSpeciesPhaseAndR2Decay(t,freqs,r2s)
%
% Input
% --------------
% t             : echo times (1xm)
% freqs         : specie frequency (1xn)
% r2s           : specie R2* (1xn)       
%
% Output
% --------------
% B             : mxn matrix with the phase and signal decay across time
%
% Description: Formula matches Nam et al. NeuroImage 2015
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 2 March 2018
% Date last modified:
%
%
function B = mwi_GetSpeciesPhaseAndR2sDecay(t,freqs,r2s)

% separated species phase and decay
B = exp(t(:)*(-r2s(:)-1i*2*pi*freqs(:)).');

end