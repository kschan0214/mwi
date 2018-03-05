%% s = mwi_model_3cm_jointT1T2s_X3lin(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1ax,t1ex,fmy,fax,kmyax,kmyex,b1,isT2sExchange)
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
function s = mwi_model_3cm_VFAT2s_2T1X_notfinish(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1ax,t1ex,fmy,fax,kmyfree,b1)

if nargin < 17
    b1=1;
end

nfa = length(fa);
nt = length(te);

kaxmy = kmyax * Amy/Aax;
kexmy = kmyex * Amy/Aex;

m0 = [Amy;Aax;Aex];
I1 = [1,1,1];

L1 = -[1/t1my,0,0;0,1/t1ax,0;0,0,1/t1ex] + [-(kmyax+kmyex),kaxmy,kexmy;kmyax,-(kaxmy+kaxex),kexax;kmyex,kaxex,-(kexmy+kexax)];
% assuming T2s << exchange time
% L2s = -[1/t2smy,1/t2sax,1/t2sex] + 1i*2*pi*[fmy,fax,0];

s = zeros(nfa,nt);
for kfa=1:nfa
    for kt=1:nt
        t1w = sind(fa(kfa)*b1) .* (eye(3)-(expm(L1 .* tr)))/(eye(3)-cosd(fa(kfa)*b1).*expm(L1.*tr)) * m0;
%         t2sw = diag(exp(L2s .* te(kt) ));
        t2sw = expm(L2s.*te(kt));

        s(kfa,kt) = abs(I1 * (t2sw * t1w));
    end
end

end