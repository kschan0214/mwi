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
% function s = mwi_model_ssSPGR_3T2scc_HCFM_MWFreq(te,A_mw,A_iw,A_ew,...
%                                                     t2smw,t2siw,t2sew,...
%                                                     freq_mw,param)
function s = mwi_model_ssSPGR_3T2scc_HCFM(te,A_mw,A_iw,A_ew,...
                                             t2smw,t2fw,...
                                             sin2theta,param)
gyro = 42.57747892;
                                     
param = checkSetDefault(param);

b0          = param.b0;
E           = param.E;
x_i         = param.x_i;
x_a         = param.x_a;
rho_mw      = param.rho_mw;
freq_bkg    = param.freq_bkg;
pini        = param.pini;
g           = param.g;
fvf         = param.fvf;

if isempty(g);      g   = sqrt(A_iw/(A_iw+A_mw/rho_mw)); end
if isempty(fvf);    fvf = (A_iw/A_ew)/((A_iw/A_ew)+g.^2); end
% s_iw = A_iw/0.91;
% s_ew = A_ew/0.8;
% s_mw = A_mw/0.37;
% g = sqrt(s_iw/(s_iw+s_mw))
% fvf = (s_iw/s_ew)./((s_iw/s_ew)+g.^2)

% coefficients, Eq.[A15]
c1 = 1/4 - (3/2)*((g^2)/(1-g^2))*log(1/g);

%
% sin2theta = sin(theta).^2;
% analytical sine square theta
% sin2theta = ((2*freq_mw/(gyro*b0)-((2*x_i-x_a)/3+2*E))/(x_a*c1-x_i));

freq_mw = ((x_i/2)*(2/3-sin2theta) + (x_a/2)*(c1*sin2theta-1/3) + E) * gyro * b0;

% analytical iw frequency
freq_iw = (3/4)*x_a*log(1/g)*sin2theta*b0*gyro;

x_d = (x_i+x_a/4)*(1-g.^2);

alpha = 3 / (abs(x_d)*2*pi*gyro*b0*sin2theta);

D_E = zeros(size(te));
D_E(te<=alpha) = (fvf/16) * (abs(x_d)*gyro*2*pi*b0*sin2theta*te(te<=alpha)).^2;
D_E(te>alpha) = (fvf/2) * (abs(x_d)*gyro*2*pi*b0*sin2theta)*(te(te>alpha)-2/((abs(x_d)*gyro*2*pi*b0*sin2theta)));

s = (A_mw*exp(te*(-1/t2smw+1i*2*pi*(freq_mw+freq_bkg))) + ...
     A_iw*exp(te*(-1/t2fw+1i*2*pi*(freq_iw+freq_bkg))) + ...
     A_ew*exp(te*(-1/t2fw+1i*2*pi*freq_bkg)).*exp(-D_E))*exp(-1i*pini);

 
end

function param2 = checkSetDefault(param)
param2 = param;

try param2.b0 = param.b0;               catch; param2.b0 = 3;       end	% T
try param2.E = param.E;                 catch; param2.E = 0.0;        end % ppm
try param2.x_i = param.x_i;             catch; param2.x_i = -0.1;   end % ppm
try param2.x_a = param.x_a;             catch; param2.x_a = -0.1;   end % ppm
try param2.rho_mw = param.rho_mw;       catch; param2.rho_mw = 0.43; end % ratio
try param2.freq_bkg = param.freq_bkg;	catch; param2.freq_bkg = 0; end % Hz
try param2.pini = param.pini;           catch; param2.pini = 0;     end % rad
try param2.g = param.g;                 catch; param2.g = [];     end % ratio
try param2.fvf = param.fvf;             catch; param2.fvf = [];     end % ratio

end