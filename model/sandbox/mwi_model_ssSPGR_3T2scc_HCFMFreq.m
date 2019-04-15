%% s = mwi_model_ssSPGR_3T2scc_HCFM(te,A_mw,A_iw,A_ew,t2smw,t2fw,sin2theta,param)
%
% Input (initial guess,[lb,ub])
% --------------
% te            : echo times
% Amy           : Myelin water signal;
% Aax           : Axonal water signal;
% Aex           : extracellular water signal;
% t2smy         : myelin water T2*;
% t2fw          : free water T2, assuming intra-axonal and extracellular water have the same intrinsic T2;
% sin2theta     : sin(theta)^2, where theta is the angle between directions
%                 of fibre and B0 field
% param         : optional parameters
%
% Output
% --------------
% s             : complex-valued 3-pool signal
%
% Description: Complex-valued model fitted to complex-valued data used in
% Nam et al. NeuroImages 2015 116:214-221
% Protocol (3T):
%   - voxel size    = 2mm isotropic
%   - fa            = 30
%   - TR            = 120ms
%   - TE            = 2.1:1.93:61.93ms (32 echoes)
%   - BW            = 1502 Hz/pixel
% Key processing:
%   (1) averaged two adjacent slices, weighted ls fitting with weight of
%   magn (voxel-by-voxel)
%   (2) averaged magn., angle of ROI-averaged complex data (ROI)
% Results (perp./para.) using 16 echoes:
%   - MWF           = 14.4/8.5 %
%   - f_(my-ex)     = 11.1/2.5 Hz
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 26 Feb 2019
% Date last modified: 
%
function s = mwi_model_ssSPGR_3T2scc_HCFMFreq(te,A_mw,A_iw,A_ew,t2smw,t2siw,t2sew,sin2theta,param)

% gyromagnetic ratio
gyro = 42.57747892;

if nargin < 8
    param = [];
end

% check and set default constants
param = checkSetDefault(param);
b0          = param.b0;
E           = param.E;
x_i         = param.x_i;
x_a         = param.x_a;
rho_mw      = param.rho_mw;
freq_bkg    = param.freq_bkg;
pini        = param.pini;
g           = param.g;
freq_mw     = param.freq_mw;

% adjust volume fractions based on water proton density
% A_iw = A_iw/0.91;
% A_ew = A_ew/0.8;
% A_mw = A_mw/0.37;
if isempty(g)
    g   = sqrt(A_iw/(A_iw+A_mw/rho_mw));
end

% coefficients, Eq.[A15]
c1 = 1/4 - (3/2)*((g^2)/(1-g^2))*log(1/g);

% analytical sin(theta)^2
if ~isempty(freq_mw)
    sin2theta = ((2*freq_mw/(gyro*b0)-((2*x_i-x_a)/3+2*E))/(x_a*c1-x_i));
else
    % analytical myelin water frequency
    freq_mw = ((x_i/2)*(2/3-sin2theta) + (x_a/2)*(c1*sin2theta-1/3) + E) * gyro * b0;
end

% analytical intra-axonal frequency
freq_iw = (3/4)*x_a*log(1/g)*sin2theta*b0*gyro;

s = (A_mw*exp(te*(-1/t2smw +1i*2*pi*(freq_mw+freq_bkg))) + ...
     A_iw*exp(te*(-1/t2siw +1i*2*pi*(freq_iw+freq_bkg))) + ...
     A_ew*exp(te*(-1/t2sew +1i*2*pi*freq_bkg)))*exp(-1i*pini);
 
end

function param2 = checkSetDefault(param)
param2 = param;

try param2.b0       = param.b0;         catch; param2.b0 = 3;      	end	% T
try param2.E        = param.E;         	catch; param2.E = 0.0;    	end % ppm
try param2.x_i      = param.x_i;      	catch; param2.x_i = -0.1;  	end % ppm
try param2.x_a      = param.x_a;     	catch; param2.x_a = -0.1;  	end % ppm
try param2.rho_mw   = param.rho_mw;  	catch; param2.rho_mw = 0.43;end % ratio
try param2.freq_bkg = param.freq_bkg;	catch; param2.freq_bkg = 0;	end % Hz
try param2.pini     = param.pini;     	catch; param2.pini = 0;   	end % rad
try param2.g        = param.g;       	catch; param2.g = [];    	end % ratio
try param2.freq_mw 	= param.freq_mw;  	catch; param2.freq_mw = []; end % Hz

end