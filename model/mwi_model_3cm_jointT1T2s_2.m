%% s = mwi_model_3cm_jointT1T2s(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1ax,t1ex,fmy,fax,b1)
%
% Input
% --------------
% fa            : flip angles
% te            : echo times
% tr            : repetition time
% Amy           : Myelin water signal
% Aax           : Axonal water signal
% Aex           : extracellular water signal
% t2smy         : myelin water T2*
% t2sax         : axonal water T2*
% t2sex         : extracellular water T2*
% t1my          : myelin water T1
% t1ax          : axonal water T1
% t1ex          : extracellular water T1
% fmy           : myelin water frequency - extracellular water
% fax           : axonal water frequency - extracellular water
% b1            : factor of B1+ inhomogeneity: true flip angle=1
%
% Output
% --------------
% s             : 2D magnitude 3-pool signal, 1st dim=fa;2nd dim=te
%
% Description: extended magnitude signal model account for the effects
% of T1w
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 18 January 2018
% Date last modified:
%
%
function s = mwi_model_3cm_jointT1T2s_2(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1l,fmy,fax,b1)

if nargin < 13
    b1=1;
end

[fa,te] = ndgrid(fa,te);

% s = sum(s_n) = sum(rho_n * T1w_n * T2s_n);
s = abs(Amy*(sind(fa*b1)*(1-exp(-tr/t1my))./(1-cosd(fa*b1)*exp(-tr/t1my))).*exp(-te*(1/t2smy+1i*2*pi*fmy)) + ...
        Aax*(sind(fa*b1)*(1-exp(-tr/t1l))./(1-cosd(fa*b1)*exp(-tr/t1l))).*exp(-te*(1/t2sax+1i*2*pi*fax)) + ...
        Aex*(sind(fa*b1)*(1-exp(-tr/t1l))./(1-cosd(fa*b1)*exp(-tr/t1l))).*exp(-te*(1/t2sex)));

end