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
function s = mwi_model_3cm_jointT1T2s_MT_test(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1ax,t1ex,fmy,fax,kmyax,kmyex,b1)

if nargin < 17
    b1=1;
end

nfa = length(fa);
nt = length(te);

kaxmy = kmyax * Amy/Aax;
kexmy = kmyex * Amy/Aex;
kaxex = 0;
kexax = 0;
% kaxex = 100;
% kexax = kaxex*Aax/Aex;

m0 = [Amy;Aax;Aex];
I1 = [1,1,1];
% Fs = 1i*2*pi*[fmy,0,0;0,fax,0;0,0,0];

L1 = -[1/t1my,0,0;0,1/t1ax,0;0,0,1/t1ex] + [-(kmyax+kmyex),kaxmy,kexmy;kmyax,-(kaxmy+kaxex),kexax;kmyex,kaxex,-(kexmy+kexax)];
maskL1 = L1~=0;
% L2s = -[1/t2smy,0,0;0,1/t2sax,0;0,0,1/t2sex] + [-(kmyax+kmyex),kaxmy,kexmy;kmyax,-(kaxmy+kaxex),kexax;kmyex,kaxex,-(kexmy+kexax)] + Fs;
% assuming T2s << exchange time
L2s = -[1/t2smy,1/t2sax,1/t2sex] + 1i*2*pi*[fmy,fax,0];

s = zeros(nfa,nt);
for kfa=1:nfa
    for kt=1:nt
        t1w = sind(fa(kfa)*b1) .* (1-(exp(L1 .* tr).*maskL1))/(1-cosd(fa(kfa)*b1).*(exp(L1.*tr).*maskL1)) * m0;
        t2sw = diag(exp(L2s .* te(kt) ));

        s(kfa,kt) = abs(I1 * (t2sw * t1w));
    end
end

% s = sum(s_n) = sum(rho_n * T1w_n * T2s_n);
% s = abs(Amy*(sind(fa*b1)*(1-exp(-tr/t1my))./(1-cosd(fa*b1)*exp(-tr/t1my))).*exp(-te*(1/t2smy+1i*2*pi*fmy)) + ...
%         Aax*(sind(fa*b1)*(1-exp(-tr/t1l))./(1-cosd(fa*b1)*exp(-tr/t1l))).*exp(-te*(1/t2sax+1i*2*pi*fax)) + ...
%         Aex*(sind(fa*b1)*(1-exp(-tr/t1l))./(1-cosd(fa*b1)*exp(-tr/t1l))).*exp(-te*(1/t2sex)));

end