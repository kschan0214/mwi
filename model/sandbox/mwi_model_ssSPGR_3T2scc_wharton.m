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
function s = mwi_model_ssSPGR_3T2scc_wharton(te,s0,fvf,g,t2smw,t2siw,t2sew,...
                                                         fmwbg,fiwbg,fbg,...
                                                         pini,rho_mw,rho_iw,rho_ew)

if nargin < 14
    rho_ew = 0.86; % Jung et al 2017 Neuroimage
end
if nargin < 13
    rho_iw = 0.86;
end
if nargin < 12
    rho_mw = 0.36;
end

Vmw = fvf*(1-g.^2);
Viw = fvf*(g.^2);
Vew = 1-fvf;

s = s0*(Vmw*rho_mw*exp(te*(-1/t2smw+1i*2*pi*fmwbg)) + ...
        Viw*rho_iw*exp(te*(-1/t2siw+1i*2*pi*fiwbg)) + ...
        Vew*rho_ew*exp(te*(-1/t2sew+1i*2*pi*fbg)))*exp(-1i*pini);

end