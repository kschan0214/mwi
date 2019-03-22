%% [freq_mw,freq_iw,theta] = FreqShift_HCFM_GivenV1(v1,param)
%
% Input
% --------------
% v1            : fibre oreintation map (cartisian coordinate system)
% param.
% ------
% x_i           : isotropic susceptibility
% x_a           : anisotropic susceptibility
% E             : chemical exchange
% g             : fibre g-ratio (ro/ri)
% b0            : field strength
% b0dir         : main field direction
%
% Output
% --------------
% freq_mw       : frequency shift of myelin water (Hz)
% freq_iw       : frequency shift of intra-axonal water (Hz)
% theta         : angle between b0 direction and fibre orientation
%
% Description: prediction of freuqnecy shift of myelin water using hollow
% cylinder model
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 26 February 2018
% Date last modified: 09 March 2018
% Date last modified: 20 March 2019
%
%
function [freq_mw,freq_iw,theta] = HCFM_FreqShift_GivenV1(v1,param)

if nargin < 2
    param = [];
end
% check and set default constants
param = checkSetDefault(param);
b0          = param.b0;
b0dir       = param.b0dir;
E           = param.E;
x_i         = param.x_i;
x_a         = param.x_a;
g           = param.g;

% compute theta
theta = AngleBetweenV1MapAndB0(v1,b0dir);

% compute frequency shift
[freq_mw,freq_iw] = HCFM_FreqShift(theta,x_i,x_a,g,E,b0);

end

function param2 = checkSetDefault(param)
param2 = param;

try param2.b0       = param.b0;         catch; param2.b0 = 3;           end	% T
try param2.b0dir    = param.b0dir;      catch; param2.b0dir = [0,0,1];  end % [x,y,z]
try param2.E        = param.E;         	catch; param2.E = 0.02;         end % ppm
try param2.x_i      = param.x_i;      	catch; param2.x_i = -0.1;       end % ppm
try param2.x_a      = param.x_a;     	catch; param2.x_a = -0.1;       end % ppm
try param2.g        = param.g;       	catch; param2.g = 0.8;          end % ratio

end