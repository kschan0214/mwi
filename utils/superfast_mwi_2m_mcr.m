%% [m0,mwf] = superfast_mwi_2m_mcr(img,te,fa,tr,t2s,t1,mask,b1map)
%
% Input
% --------------
% img           : variable flip angle multi-echo GRE image, 5D [row,col,slice,TE,Flip angle]
% te            : echo times in second
% fa            : flip angle in degree
% tr            : repetition time in second
% t2s           : T2* of the two pools, in second, [T2sMW,T2sIEW], if empty
%                 then default values for 3T will be used
% t1            : T1 of the two pools, in second, [T1MW, T1IEW], if empty
%                 then default values for 3T will be used
% mask          : signal mask, (optional)
% b1map         : B1 flip angel ratio map, (optional)
%
% Output
% --------------
% m0            : proton density of each pool, 4D [row,col,slice,pool]
% mwf           : myelin water fraction map, range [0,1]
%
% Description:  Direct matrix inversion based on simple 2-pool model, i.e.
%               S(te,fa) = E1 * M0 * E2s
%               Useful to estimate initial starting points for MWI fitting
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 Nov 2020
% Date modified:
%
%
function [m0,mwf] = superfast_mwi_2m_mcr(img,te,fa,tr,t2s,t1,mask,b1map)

% get size in all image dimensions
dims(1) = size(img,1);
dims(2) = size(img,2);
dims(3) = size(img,3);

% check input
if isempty(t2s)
    t2s = [10e-3, 60e-3];   % 3T, [MW, IEW], in second
end
if isempty(t1)
    t1 = [234e-3, 1];       % 3T, [MW, IEW], in second
end
if nargin < 8
    b1map = ones(dims);
end
if nargin < 7
    mask = ones(dims);
end

% T2* decay matrix
E2s1    = exp(-te(:).'/t2s(1));
E2s2	= exp(-te(:).'/t2s(2));
E2s     = [E2s1;E2s2];

tmp = reshape(abs(img),prod(dims),length(te),length(fa));

m0   = zeros([prod(dims),2]);
for k = 1:prod(dims)
    if mask(k) ~= 0
        
        % T1-T2* signal
        temp = squeeze(tmp(k,:,:)).';

        % T1 steady-state matrix
        E11 = sind(fa(:)*b1map(k)) .* (1-exp(-tr/t1(1)))./(1-cosd(fa(:)*b1map(k))*exp(-tr/t1(1)));
        E12 = sind(fa(:)*b1map(k)) .* (1-exp(-tr/t1(2)))./(1-cosd(fa(:)*b1map(k))*exp(-tr/t1(2)));
        E1  = [E11,E12];
        
        % matrix inversion
        temp2 = (E1 \ temp) / E2s;
        
        % the diagonal components represent signal amplitude
        m0(k,:) =  diag(temp2);
    end

end
m0 = reshape(m0,[dims 2]);

% compute MWF
mwf = m0(:,:,:,1) ./ sum(m0,4);
mwf(mwf<0)      = 0;
mwf(isnan(mwf)) = 0;
mwf(isinf(mwf)) = 0;

m0(m0 < 0)      = 0;
m0(isinf(m0))   = 0;
m0(isnan(m0))   = 0;

end