%% function output = function_name(input)
%
% Input
% --------------
% x_i           : isotropic susceptibility (in ppm)
% x_a           : anisotropic susceptibility (in ppm)
% g             : g-ratio
% theta         : angle between B0 and fibre direction (degree)
% E             : exchange (in ppm)
% b0            : field strength (in T)
%
% Output
% --------------
% freq_myelin   : frequency shift in myelin (in Hz)
% freq_axon     : frequency shift in axon (in Hz)
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 July 2019
% Date last modified:
%
%
function [freq_myelin, freq_axon] = hcfm_freq(x_i,x_a,g,theta,E,b0)

gyro = 42.57747892;

if ~exist('b0','var') || isempty(b0);  b0  = 3; end
if ~exist('E','var') || isempty(E);    E   = 0.02; end

% Eq. [A15]
c1 = hcfm_c1(g);

% Eq.[5]
freq_myelin = ((x_i./2).*(2/3 - sind(theta).^2) + (x_a./2).*(c1.*sind(theta).^2 - 1/3) + E)*gyro*b0;
% Eq.[6]
freq_axon = (3*x_a./4) .* sind(theta).^2 .* log(1./g) * gyro*b0;

end