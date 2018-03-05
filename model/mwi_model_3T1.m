%% s = mwi_model_3T1(fa,tr,Amy,Aax,Aex,t1my,t1ax,t1ex,b1)
%
% Input (initial guess,[lb,ub])
% --------------
% fa            : flip angles
% tr            : repetition time
% Amy           : Myelin water signal;          (0.1*abs(S1),[0,2*abs(S1)])
% Aax           : Axonal water signal;         	(0.6*abs(S1),[0,2*abs(S1)])
% Aex           : extracellular water signal;	(0.3*abs(S1),[0,2*abs(S1)])
% t1my          : myelin water T1;            	(537,[100,650])
% t1ax          : axonal water T1;             	(1183,[650,2000])
% t1ex          : extracellular water T1;      	(1183,[650,2000])
% b1            : B1+ ratio
%
% Output
% --------------
% s             : magnitude-valued 3-pool signal
%
% Description: 
% Key processing:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 2 March 2018
% Date last modified:
%
%
function s = mwi_model_3T1(fa,tr,Amy,Aax,Aex,t1my,t1ax,t1ex,b1)

if nargin < 9
    b1=1;
end

fa = b1*fa;

E1my = exp(-tr./t1my);
E1ax = exp(-tr./t1ax);
E1ex = exp(-tr./t1ex);
s = (Amy.*(1-E1my).*sind(fa))./(1-E1my.*cosd(fa)) + ...
    (Aax.*(1-E1ax).*sind(fa))./(1-E1ax.*cosd(fa)) + ...
    (Aex.*(1-E1ex).*sind(fa))./(1-E1ex.*cosd(fa));

end