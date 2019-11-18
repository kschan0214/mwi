%% s = mwi_model_3T1_ssSPGR(fa,tr,m0mw,m0iw,m0ew,t1mw,t1iw,t1ew)
%
% Input
% --------------
% fa            : flip angle, in degree
% tr            : repetition time, in second
% m0mw      	: myelin water proton density
% m0iw       	: intra-aonxal water proton density
% m0ew          : extracellular water proton density
% t1mw          : myelin water T1, in second
% t1iw          : intra-aonxal water T1, in second
% t1ew          : extracellular water T1, in second
%
% Output
% --------------
% s             :  3-pool T1w signal
%
% Description: 
% Date created: 15 November 2019
% Date last modified: 
%
%
function ss = mwi_model_3T1_ssSPGR(fa,tr,m0mw,m0iw,m0ew,t1mw,t1iw,t1ew,b1)

if nargin < 9
    b1 = 1;
end

ss = zeros(3,length(fa));

fa = fa * b1;

ss(1,:) = m0mw .* sind(fa) .* (1-exp(-tr/t1mw))./(1-cosd(fa).*exp(-tr/t1mw));
ss(2,:) = m0iw .* sind(fa) .* (1-exp(-tr/t1iw))./(1-cosd(fa).*exp(-tr/t1iw));
ss(3,:) = m0ew .* sind(fa) .* (1-exp(-tr/t1ew))./(1-cosd(fa).*exp(-tr/t1ew));
    
    
end