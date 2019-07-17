%% [v_myelin,v_axon,v_ec] = hcfm_volume_fraction(fvf,g)
%
% Input
% --------------
% fvf           : fibre volume fraction
% g             : g-ratio
%
% Output
% --------------
% v_myelin      : volume fraction of myelin
% v_axon        : volume fraction of axon
% v_ec          : volume fraction of extracellular space
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 July 2019
% Date last modified:
%
%
function [v_myelin,v_axon,v_ec] = hcfm_volume_fraction(fvf,g)

v_myelin	= fvf .* (1-g.^2);
v_axon  	= fvf .* g.^2;
v_ec        = 1-fvf;

end