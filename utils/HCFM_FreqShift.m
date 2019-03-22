%% [freq_mw,freq_iw] = HCFM_FreqShift(theta,x_i,x_a,g,E,b0)
%
% Input
% --------------
% theta         : angle between fibre orientation and b0 direction (rad)
% x_i           : isotropic susceptibility of myelin (ppm)
% x_a           : anisotropic susceptibility of myelin (ppm)
% g             : g-ratio 
% E             : exchange term (ppm)
% b0            : magnetic field strength (T)
%
% Output
% --------------
% freq_mw       : frequency shift in myelin shealth
% freq_iw       : frequency shift in intra-axonal space
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 20 March 2019
% Date last modified:
%
%
function [freq_mw,freq_iw] = HCFM_FreqShift(theta,x_i,x_a,g,E,b0)
% gyromagnetic ratio (MHz/T)
gyro = 42.57747892 ;

sin2theta = sin(theta).^2;
c1 = 1/4 - (3/2)*((g.^2)/(1-g.^2))*log(1/g);

freq_mw = ((x_i/2)*(2/3-sin2theta) + (x_a/2)*(c1*sin2theta-1/3) + E) * gyro * b0;
freq_iw = (3/4)*x_a*log(1/g)*sin2theta*gyro*b0;

end