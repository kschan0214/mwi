%% fiter = computeFiter(s,shat,NUM_MAGN,w)
% s     - measured signal
% shat  - simulated signal
% NUM_GAGN - no. of phase corrupted echoes:
% NUM_MAGN=0 : complex fitting
% NUM_MAGN=length(s) : magnitude fitting
% NUM_MAGN (0,length(s)) : mixed fitting
% w : weights, must be same size as s
%
% Description: Compute the fitter for lsqnonlin
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
function fiter = computeFiter(s,shat,NUM_MAGN,w)
if nargin < 4
    w = ones(size(s));
end

if isvector(s)
    nt = length(s);
else
    % if dimension of s > 2, then assume time is in the second dimenison
    nt = size(s,2);
end

% compute the residual
if NUM_MAGN == nt
%% Option 1: Magnitude fitting, when the number of phase corrupted-echoes == length of s
    % compute the magnitude of both measurement and estimation
    shat1 = abs(shat);
    s1 = abs(s);
    % apply weights
    fiter = w .* (shat1-s1);
    % vectorise the residual
    fiter = fiter(:);
    
    % deprectaed
%     fiter = shat1(:) - s1(:);
%     w = w(:);

elseif NUM_MAGN == 0 
%% Option 2: Complex fitting, when the no. of phase-corrupted echoes == 0
    % compute the residual of complex-valued data and apply weights
    fiter = w .* (shat-s);
    % separate real and imaginary parts and vectorise the residual
    fiter = [real(fiter(:)); imag(fiter(:))];
    
    % deprecated
%         fiter2 = shat(:) - s(:);
%         fiter2 = [real(fiter2); imag(fiter2)];
%         fiter2 = [real(fiter2), imag(fiter2)];
%         fiter = fiter2;
%         w = repmat(w(:),2,1);

else
%% Option 3: Mixed fitting
    % 1D measurement
    if isvector(s)
        w = w(:);
        
        shat1 = abs(shat(1:NUM_MAGN));
        s1 = abs(s(1:NUM_MAGN));
        shat2 = shat(NUM_MAGN+1:end);
        s2 = s(NUM_MAGN+1:end);


        fiter1 = w(1:NUM_MAGN) .* (shat1(:)-s1(:));
        fiter2 = w(NUM_MAGN+1:end) .* (shat2(:)-s2(:));
        fiter2 = [real(fiter2(:));imag(fiter2(:))];

        fiter = [fiter1(:);fiter2];
        
        % deprecated
%         fiter1 = shat1(:) - s1(:);
%         fiter2 = shat2(:) - s2(:);
%         fiter2 = [real(fiter2); imag(fiter2)];
%         fiter = [fiter1(:);fiter2];        
%         w1 = w(1:NUM_MAGN);
%         w2 = w(NUM_MAGN+1:end);
%         w = [w1(:);repmat(w2(:),2,1)];
    else
        % 2D measurement
        shat1 = abs(shat(:,1:NUM_MAGN));
        s1 = abs(s(:,1:NUM_MAGN));
        shat2 = shat(:,NUM_MAGN+1:end);
        s2 = s(:,NUM_MAGN+1:end);

        fiter1 = w(:,1:NUM_MAGN) .* (shat1-s1);
        fiter1 = fiter1(:);
        fiter2 = w(:,NUM_MAGN+1:end) .* (shat2-s2);
        fiter2 = [real(fiter2(:));imag(fiter2(:))];

        fiter = [fiter1;fiter2];
        
        % deprecated
%         fiter1 = shat1(:) - s1(:);
%         fiter2 = shat2(:) - s2(:);
%         fiter2 = [real(fiter2); imag(fiter2)];
%         fiter = [fiter1;fiter2];
%         w1 = w(:,1:NUM_MAGN);
%         w2 = w(:,NUM_MAGN+1:end);
%         w = [w1(:);repmat(w2(:),2,1)];
    end
    
end
fiter = double(fiter);
% fiter = double(fiter.*w);
end