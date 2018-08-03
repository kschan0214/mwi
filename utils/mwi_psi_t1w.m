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
function psi = mwi_psi_t1w(t1_vec, fa, tr)

r1_vec = 1 ./ t1_vec(:);
fa = fa(:).';

psi = ((1-exp(-tr*r1_vec)) * sind(fa)) ./ (1 - exp(-tr*r1_vec) * cosd(fa));

end