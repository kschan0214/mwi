%% freqMye = FreqShiftMyelin_HollowCylinder(v1,chiI,chiA,E,gratio,b0)
%
% Input
% --------------
% v1            : fibre oreintation map (cartisian coordinate system)
% chiI          : isotropic susceptibility
% chiA          : anisotropic susceptibility
% E             : chemical exchange
% gratio        : fibre g-ratio (ro/ri)
% b0            : field strength
%
% Output
% --------------
% freqMye       : frequency shift of myelin water
%
% Description: prediction of freuqnecy shift of myelin water using hollow
% cylinder model
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 26 February 2018
% Date last modified:
%
%
function freqMye = FreqShiftMyelin_HollowCylinder(v1,chiI,chiA,E,gratio,b0)
gyro = 42.57747892 * 1e6;

if nargin < 2
    % parameters uded in Nam's paper
    b0=3;
    chiI = -100e-9; % -100 ppb
    chiA = -100e-9; % -100 ppb
    E = 20e-9; % 20 ppb
    gratio=0.8;
end

cthe = v1(:,:,:,3)./(sqrt(sum(v1.^2,4)));
sthe = (sqrt(sum(v1(:,:,:,1:2).^2,4)))./(sqrt(sum(v1.^2,4)));

freqMye = gyro*b0 * (chiI/2*(cthe.^2-1/3) + E + chiA/2*(-1/3+(sthe.^2)*(1/4-(3/2)*(1^2/(gratio.^2-1.^2))*log(gratio))));
end