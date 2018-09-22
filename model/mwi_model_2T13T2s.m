%% s = mwi_model_2T13T2s(afa,te,tr,s0,mwf,iaf,t2smw,t2siw,t2sew,t1mw,t1ow,fmw,fiw,b1)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
%
function s = mwi_model_2T13T2s(afa,te,tr,s0,mwf,iaf,t2smw,t2siw,t2sew,t1mw,t1ow,fmw,fiw,b1)

if nargin < 15
    b1=1;
end

[te,cfa] = ndgrid(te,afa*b1);

% s = sum(s_n) = sum(rho_n * T1w_n * T2s_n);
s = s0*(mwf*            (sind(cfa)*(1-exp(-tr/t1mw))./(1-cosd(cfa)*exp(-tr/t1mw))).*exp(-te*(1/t2smw-1i*2*pi*fmw)) + ...
        (1-mwf)*iaf*    (sind(cfa)*(1-exp(-tr/t1ow))./(1-cosd(cfa)*exp(-tr/t1ow))).*exp(-te*(1/t2siw-1i*2*pi*fiw)) + ...
        (1-mwf)*(1-iaf)*(sind(cfa)*(1-exp(-tr/t1ow))./(1-cosd(cfa)*exp(-tr/t1ow))).*exp(-te*(1/t2sew)));
    
end