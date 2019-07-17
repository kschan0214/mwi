%% x_d = hcfm_effective_susceptibility(x_i,x_a,g)
%
% Input
% --------------
% x_i           : isotropic susceptibility (in ppm)
% x_a           : anisotropic susceptibility (in ppm)
% g             : g-ratio
% 
% Output
% --------------
% x_d           : effective susceptibility (in ppm)
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 July 2019
% Date last modified:
%
%
function x_d = hcfm_effective_susceptibility(x_i,x_a,g)

% Eq. [A9]
x_d = (x_i + x_a/4).*(1-g.^2);

end