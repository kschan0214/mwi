%% s = mwi_model_3cc_nam2015(te,Amy,Aax,Aex,t2smy,t2sax,t2sex,fmy,fax,fex,pini)
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
% fmybg         : myelin water frequency + background field;         (fbkg,[fbkg-75,fbkg+75]Hz)
% faxbg         : axonal water frequency + background field;         (fbkg,[fbkg-25,fbkg+25]Hz)
% fexbg         : exctracellular water frequency + background field; (fbkg,[fbkg-25,fbkg+25]Hz)
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
% Date created: 4 January 2018
% Date last modified: 12 February 2018
%
%
function s = mwi_model_3cc_nam2015(te,Amy,Aax,Aex,t2smy,t2sax,t2sex,fmybg,faxbg,fexbg,pini)

if nargin < 10
    pini=0;
end

s = (Amy*exp(te*(-1/t2smy+1i*2*pi*fmybg)) + ...
     Aax*exp(te*(-1/t2sax+1i*2*pi*faxbg)) + ...
     Aex*exp(te*(-1/t2sex+1i*2*pi*fexbg)))*exp(-1i*pini);

end