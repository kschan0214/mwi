%% D_E = hcfm_dephasing_decay_DE(t,fvf,g,x_i,x_a,theta,b0)
%
% Input
% --------------
% t             : time (in second)
% fvf           : fibre volume fraction
% g             : g-ratio
% x_i           : isotropic susceptibility (in ppm)
% x_a           : anisotropic susceptibility (in ppm)
% theta         : angle between B0 and fibre direction (rad)
% b0            : field strength (in T)
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 July 2019
% Date last modified:
%
%
function D_E = hcfm_dephasing_decay_DE(t,fvf,g,x_i,x_a,theta,b0)
gyro = 2*pi*42.57747892;

if ~exist('b0','var') || isempty(b0);  b0  = 3; end

alpha = hcfm_transition_time(x_i,x_a,g,theta,b0);

x_d = hcfm_effective_susceptibility(x_i,x_a,g);

D_E = zeros(size(t));
for kt=1:length(t)
    if t(kt) < alpha
        D_E(kt) = (fvf/16)*(abs(x_d)*gyro*b0*sin(theta).^2*t(kt))^2;
    else
        D_E(kt) = (fvf/2)*(abs(x_d)*gyro*b0*sin(theta).^2)*...
            (t(kt) - 2./(abs(x_d)*gyro*b0*sin(theta).^2));
    end
end


end