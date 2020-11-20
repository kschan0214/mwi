%% [m0,mwf] = superfast_mwi_3m_standard(img,te,t2s)
%
% Input
% --------------
% img           : variable flip angle multi-echo GRE image, 5D [row,col,slice,TE,Flip angle]
% te            : echo times in second
% t2s           : T2* of the two pools, in second, [T2sMW,T2sEW, T2sIW], if empty
%                 then default values for 3T will be used
%
% Output
% --------------
% m0            : proton density of each pool, 4D [row,col,slice,pool]
% mwf           : myelin water fraction map, range [0,1]
%
% Description:  Direct matrix inversion based on simple 2-pool model, i.e.
%               S(te,fa) = E2s * M0
%               Useful to estimate initial starting points for MWI fitting
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 Nov 2020
% Date modified:
%
%
function [m0,mwf] = superfast_mwi_3m_standard(img,te,t2s)

% get size in all image dimensions
dims(1) = size(img,1);
dims(2) = size(img,2);
dims(3) = size(img,3);

% check input
if isempty(t2s)
    t2s = [10e-3, 48e-3, 64e-3];   % 3T, [MW, EW, IW], in second
end

% T2* decay matrix
E2s1    = exp(-te(:)/t2s(1));
E2s2	= exp(-te(:)/t2s(2));
E2s3	= exp(-te(:)/t2s(3));
E2s     = [E2s1,E2s2,E2s3];

tmp = reshape(abs(img),prod(dims),length(te));

m0 = E2s \ tmp.';
m0 = reshape(m0.',[dims length(t2s)]);

% compute MWF
mwf = m0(:,:,:,1) ./ sum(m0,4);
mwf(mwf<0)      = 0;
mwf(isnan(mwf)) = 0;
mwf(isinf(mwf)) = 0;

end