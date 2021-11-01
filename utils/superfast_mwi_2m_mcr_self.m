%% [m0,mwf,t2s_iew,t1_iew] = superfast_mwi_2m_mcr_self(img,te,fa,tr,t2s,t1_mw,mask,b1map,mode)
%
% Input
% --------------
% img           : variable flip angle multi-echo GRE image, 5D [row,col,slice,TE,Flip angle]
% te            : echo times in second
% fa            : flip angle in degree
% tr            : repetition time in second
% t2s           : T2* of the two pools, in second, [T2sMW,T2sIEW], if empty
%                 then default values for 3T will be used
% t1_mw         : T1 of MW, in second, if empty
%                 then a default value for 3T will be used
% mask          : signal mask, (optional)
% b1map         : B1 flip angel ratio map, (optional)
% mode          : IEW T1 estimation approach, ('superfast' or 'normal')
%
% Output
% --------------
% m0            : proton density of each pool, 4D [row,col,slice,pool]
% mwf           : myelin water fraction map, range [0,1]
% t2s_iew       : IEW T2*, in second
% t1_iew        : IEW T1, in second
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
function [m0,mwf,t2s_iew,t1_iew] = superfast_mwi_2m_mcr_self(img,te,fa,tr,t2s,t1_mw,mask,b1map,mode)

% get size in all image dimensions
dims(1) = size(img,1);
dims(2) = size(img,2);
dims(3) = size(img,3);

% check input
if isempty(t2s)
    t2s = [10e-3, 60e-3];   % 3T, [MW, IEW], in second
end
if isempty(t1_mw)
    t1_mw  = 234e-3; % 3T
end
if nargin < 9
    mode = 'superfast';
end
if nargin < 8 || isempty(b1map)
    b1map = ones(dims);
end
if nargin < 7 || isempty(mask)
    mask = ones(dims);
end

% assign maximum T2* to IEW
ind     = find(te>1.5*t2s(1));
if isempty(ind)
    t2s_iew = zeros([dims length(fa)]);
    for kk = 1:length(fa)
        [~,t2s_iew(:,:,:,kk),~] = R2star_trapezoidal(abs(img(:,:,:,:,kk)),te);
    end
else
    t2s_iew = zeros([dims length(te)-length(ind) length(fa)]);
    for kk = 1:length(fa)
        for k = 1:length(te)-length(ind)
            [~,t2s_iew(:,:,:,k,kk),~] = R2star_trapezoidal(abs(img(:,:,:,k:end,kk)),te);
        end
    end
end
t2s_iew = reshape(t2s_iew,[dims size(t2s_iew,4)*size(t2s_iew,5)]);
t2s_iew = max(t2s_iew,[],4);

t2s_iew = reshape(t2s_iew,prod(dims),1);

t2s_mw = t2s(1);

switch mode
    case 'superfast'

        % T2* decay matrix
        E2s1    = exp(-te(:)/t2s(1));
        E2s2	= exp(-te(:)/t2s(2));
        E2s     = [E2s1,E2s2];

%         s0_mw   = zeros([dims length(fa)]);
        s0_iew  = zeros([dims length(fa)]);
        for kfa = 1:length(fa)

            tmp = reshape(abs(img(:,:,:,:,kfa)),prod(dims),length(te));

            s0 = E2s \ tmp.';

%             s0_mw(:,:,:,kfa)    = reshape(s0(1,:),dims);
            s0_iew(:,:,:,kfa)   = reshape(s0(2,:),dims);
        end
        % IEW T1 estimation is relative insensitive to the change of its T2*
        [t1_iew, ~] = despot1_mapping(s0_iew,fa,tr,mask,b1map);
        
    otherwise
        
        % get IEW signal amplitude for DESPOT1, which should be more robust
        % s0_mw   = zeros([dims length(fa)]);
        s0_iew  = zeros([dims length(fa)]);
        for kfa = 1:length(fa)

        tmp = reshape(abs(img(:,:,:,:,kfa)),prod(dims),length(te));

        s0 = zeros(2,prod(dims));
        for k = 1:prod(dims)

            % T2* decay matrix
            E2s = [exp(-te(:)/t2s_mw),exp(-te(:)/t2s_iew(k))];

            s0(:,k) = E2s \ tmp(k,:).';
        end
        % s0_mw(:,:,:,kfa)    = reshape(s0(1,:),dims);
        s0_iew(:,:,:,kfa)   = reshape(s0(2,:),dims);
        end
        % compute IEW T1
        % [t1_mw, m0_mw]      = despot1_mapping(s0_mw,fa,tr,mask,b1map);
        [t1_iew, ~]    = despot1_mapping(s0_iew,fa,tr,mask,b1map);
end

% main
tmp = reshape(abs(img),prod(dims),length(te),length(fa));
m0   = zeros([prod(dims),2]);
for k = 1:prod(dims)
    if mask(k) ~= 0
    
        % T1-T2* signal
        temp = squeeze(tmp(k,:,:)).';

        % T2* decay matrix
        E2s = [exp(-te(:).'/t2s(1));exp(-te(:).'/t2s_iew(k))];

        % T1 steady-state matrix
        E1_mw  = sind(fa(:)*b1map(k)) .* (1-exp(-tr/t1_mw))./(1-cosd(fa(:)*b1map(k))*exp(-tr/t1_mw));
        E1_iew = sind(fa(:)*b1map(k)) .* (1-exp(-tr/t1_iew(k)))./(1-cosd(fa(:)*b1map(k))*exp(-tr/t1_iew(k)));
        E1 = [E1_mw,E1_iew];
        
        % matrix inversion
%         temp2 = pinv(E1)*temp*pinv(E2s);
        temp2 = (E1 \ temp) / (E2s);
        
        % the diagonal components represent signal amplitude
        m0(k,:) =  diag(temp2);
    end
end

t2s_iew = reshape(t2s_iew,dims);

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