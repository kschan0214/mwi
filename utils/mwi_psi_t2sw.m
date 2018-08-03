%% psi = mwi_psi_t2sw(t2s_vec, freq_vec, t)
%
% Input
% --------------
% t2s_vec       : 1-by-nspecie vector of T2*s, where n is the no. of species
% freq_vec      : 1-by-nspecie vector of frequncy shifts
% t             : echo times with nt echoes
%
% Output
% --------------
% psi           : nt-by-nspecie T2* weighting and frequency shifting
%                 effects
%
% Description: T2* weighting of individual species given their T2*s and
% frequency shifts at times t
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 May 2018
% Date last modified:
%
%
function psi = mwi_psi_t2sw(t2s_vec, freq_vec, t)

r2s_vec = 1 ./ t2s_vec(:).';
freq_vec = freq_vec(:).';
t = t(:);

psi = exp(t * (-r2s_vec + 1i*2*pi*freq_vec));

end