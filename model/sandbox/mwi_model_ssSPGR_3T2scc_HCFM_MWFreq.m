%% s = mwi_model_ssSPGR_3T2scc_wharton(te,s0,fvf,g,rho_mw,rho_iw,rho_ew,...
%                                                  t2smw,t2siw,t2sew,...
%                                                  fmwbg,fiwbg,fbg,pini)
%
% Input (initial guess,[lb,ub])
% --------------
% te            : echo times
% Amy           : Myelin water signal;                               (0.1*abs(S1),[0,2*abs(S1)])
% Aax           : Axonal water signal;                               (0.6*abs(S1),[0,2*abs(S1)])
% Aex           : extracellular water signal;                        (0.3*abs(S1),[0,2*abs(S1)])
% t2smy         : myelin water T2*;                                  (10,[3,25])
% t2sax         : axonal water T2*;                                  (64,[25,150])
% t2sex         : extracellular water T2*;                           (48,[25,150])
% fmwbg         : myelin water frequency + background field;         (fbkg,[fbkg-75,fbkg+75]Hz)
% fiwbg         : axonal water frequency + background field;         (fbkg,[fbkg-25,fbkg+25]Hz)
% fbg           : background field;                                  (fbkg,[fbkg-25,fbkg+25]Hz)
% pini          : initial phase introduced by B1+ phase offset;      (angle(S1),[-pi,pi])
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
%
function s = mwi_model_ssSPGR_3T2scc_HCFM_MWFreq(te,A_mw,A_iw,A_ew,...
                                         t2smw,t2siw,t2sew,...
                                         freq_mw,...
                                         param)
gyro = 42.57747892;
                                     
param = checkSetDefault(param);

b0 = param.b0;
E = param.E;
x_i = param.x_i;
x_a = param.x_a;
rho_mw = param.rho_mw;

g = sqrt(A_iw/(A_iw+A_mw/rho_mw));
% coefficients, Eq.[A15]
c1 = 1/4 - (3/2)*((g^2)/(1-g^2))*log(1/g);

% analytical iw frequency
freq_iw = (3/4)*x_a*log(1/g)*((2*freq_mw/(gyro*b0)-((2*x_i-x_a)/3+2*E))/(x_a*c1-x_i))*b0*gyro;


s = (A_mw*exp(te*(-1/t2smw+1i*2*pi*(freq_mw+freq_bkg))) + ...
     A_iw*exp(te*(-1/t2siw+1i*2*pi*(freq_iw+freq_bkg))) + ...
     A_ew*exp(te*(-1/t2sew+1i*2*pi*freq_bkg)))*exp(-1i*pini);

 
end

function param2 = checkSetDefault(param)
param2 = param;

try param2.b0 = param.b0;               catch; param2.b0 = 3;       end	% T
try param2.E = param.E;                 catch; param2.E = 0;        end % ppm
try param2.x_i = param.x_i;             catch; param2.x_i = -0.1;   end % ppm
try param2.x_a = param.x_a;             catch; param2.x_a = -0.1;   end % ppm
try param2.rho_mw = param.rho_mw;       catch; param2.rho_mw = 0.7; end % ratio
try param2.freq_bkg = param.freq_bkg;	catch; param2.freq_bkg = 0; end % Hz
try param2.pini = param.pini;           catch; param2.pini = 0;     end % rad

end