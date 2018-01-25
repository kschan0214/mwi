%% s = mwi_model_3cm_vGelderen2012(te,Amy,Aax,Aex,t2smy,t2sax,t2sex,fmy,fax,fex,pini)
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
% fmy           : myelin water frequency - extracellular water;      (11.1,[11.1-75,11.1+75]Hz)
% fax           : axonal water frequency - extracellular water;      (0,[0-25,0+25]Hz)
%
% Output
% --------------
% s             : complex-valued 3-pool signal
%
% Description: Complex-valued model fitted to magnitude data used in
% Nam et al. NeuroImages 2015 116:214-221 and van Gelderen MRM 2012
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
% Date created: 18 January 2018
% Date last modified:
%
%
function s = mwi_model_3cm_vGelderen2012(te,Amy,Aax,Aex,t2smy,t2sax,t2sex,fmy,fax)

s = abs(Amy*exp(-te*(1/t2smy+1i*2*pi*fmy)) + ...
        Aax*exp(-te*(1/t2sax+1i*2*pi*fax)) + ...
        Aex*exp(-te*(1/t2sex)));

end