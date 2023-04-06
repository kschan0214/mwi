%% [m0,mwf,t2s_iew] = superfast_mwi_2m_standard_self(img,te,t2s)
%
% Input
% --------------
% img           : variable flip angle multi-echo GRE image, 5D [row,col,slice,TE,Flip angle]
% te            : echo times in second
% t2s           : T2* of the two pools, in second, [T2sMW,T2sIEW], if empty
%                 then default values for 3T will be used
%
% Output
% --------------
% m0            : proton density of each pool, 4D [row,col,slice,pool]
% mwf           : myelin water fraction map, range [0,1]
% t2s_iew       : IEW T2*, in second
%
% Description:  Direct matrix inversion based on simple 2-pool model, i.e.
%               S(te) = E2s * M0 
%               Useful to estimate initial starting points for MWI fitting
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 Nov 2020
% Date modified:
%
%
function [m0,mwf,t2s_iew] = superfast_mwi_2m_standard_self(img,te,t2s_mw)

% get size in all image dimensions
dims(1) = size(img,1);
dims(2) = size(img,2);
dims(3) = size(img,3);

% check input
if isempty(t2s_mw)
    t2s_mw = 10e-3;   % 3T, [MW, IEW], in second
end

% assign maximum T2* to IEW
ind     = find(te>1.5*t2s_mw);
if isempty(ind)
    [~,t2s_iew,~] = R2star_trapezoidal(abs(img),te);
else
    t2s_iew = zeros([dims length(te)-length(ind)]);
    for k = 1:length(te)-length(ind)
        [~,t2s_iew(:,:,:,k),~] = R2star_trapezoidal(abs(img(:,:,:,k:end)),te);
    end
    t2s_iew = max(t2s_iew,[],4);
end
t2s_iew = reshape(t2s_iew,prod(dims),1);

tmp = reshape(abs(img),prod(dims),length(te));
m0 = zeros(2,prod(dims));
for k = 1:prod(dims)

    % T2* decay matrix
    E2s = [exp(-te(:)/t2s_mw),exp(-te(:)/t2s_iew(k))];

    m0(:,k) = E2s \ tmp(k,:).';
end

t2s_iew = reshape(t2s_iew,dims);

m0 = reshape(m0.',[dims 2]);

% compute MWF
mwf = m0(:,:,:,1) ./ sum(m0,4);
mwf(mwf<0)      = 0;
mwf(mwf>1)      = 1;
mwf(isnan(mwf)) = 0;
mwf(isinf(mwf)) = 0;

m0(m0 < 0)      = 0;
m0(isinf(m0))   = 0;
m0(isnan(m0))   = 0;

end