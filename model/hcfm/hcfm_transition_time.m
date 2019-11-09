%% t = hcfm_transition_time(x_i,x_a,g,theta,b0)
%
% Input
% --------------
% x_i           : isotropic susceptibility (in ppm)
% x_a           : anisotropic susceptibility (in ppm)
% g             : g-ratio
% theta         : angle between B0 and fibre direction (rad)
% b0            : field strength (in T)
%
% Output
% --------------
% t             : transition time of dephasing from quadratic regime to
%                 linear regime
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 July 2019
% Date last modified:
%
%
function t = hcfm_transition_time(x_i,x_a,g,theta,b0)
gyro = 2*pi*42.57747892;

if ~exist('b0','var') || isempty(b0);  b0  = 3; end

% Eq. A9
x_d = hcfm_effective_susceptibility(x_i,x_a,g);

% Eq. A8
t = 3 ./ (abs(x_d)*gyro*b0.*sin(theta).^2);

end