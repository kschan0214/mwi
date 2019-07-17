%% c1 = hcfm_c1(g)
%
% Input
% --------------
% g             : g-ratio
%
% Output
% --------------
% c1            : c1 of equation A15
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 July 2019
% Date last modified:
%
%
function c1 = hcfm_c1(g)
c1 = 1/4 - (3/2)*((g.^2)./(1-g.^2)).*log(1./(g.^2));
c1(g==1) = 1/4 - 3/2;
end