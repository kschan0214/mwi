%% function fiter = computeFiter(s,shat,NUM_MAGN)
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
    nt = size(s,2);
end

if NUM_MAGN == nt % Magnitude fitting
    shat1 = abs(shat);
    s1 = abs(s);
    fiter = shat1(:) - s1(:);
    w = w(:);
elseif NUM_MAGN == 0 % Complex fitting
        fiter2 = shat(:) - s(:);
        fiter2 = [real(fiter2); imag(fiter2)];
%         fiter2 = [real(fiter2), imag(fiter2)];
        fiter = fiter2;
        w = repmat(w(:),2,1);
else
    % Compute mixed fitting fit error
    if isvector(s)
        shat1 = abs(shat(1:NUM_MAGN));
        s1 = abs(s(1:NUM_MAGN));
        shat2 = shat(NUM_MAGN+1:end);
        s2 = s(NUM_MAGN+1:end);

        fiter1 = shat1(:) - s1(:);
        fiter2 = shat2(:) - s2(:);
        fiter2 = [real(fiter2); imag(fiter2)];

        fiter = [fiter1;fiter2];
        
        w1 = w(1:NUM_MAGN);
        w2 = w(NUM_MAGN+1:end);
        w = [w1(:);repmat(w2(:),2,1)];
    else
        shat1 = abs(shat(:,1:NUM_MAGN));
        s1 = abs(s(:,1:NUM_MAGN));
        shat2 = shat(:,NUM_MAGN+1:end);
        s2 = s(:,NUM_MAGN+1:end);

        fiter1 = shat1(:) - s1(:);
        fiter2 = shat2(:) - s2(:);
        fiter2 = [real(fiter2); imag(fiter2)];

        fiter = [fiter1;fiter2];
        
        w1 = w(:,1:NUM_MAGN);
        w2 = w(:,NUM_MAGN+1:end);
        w = [w1(:);repmat(w2(:),2,1)];
    end
    
end
fiter = double(fiter.*w);
end