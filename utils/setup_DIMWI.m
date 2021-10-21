%% [DIMWIParams, icvf, ff, theta] = setup_DIMWI(imgParam, algoParam)
%
% Input
% --------------
% algoParam     : structure array contains all fitting algorithm specific
%                 parameters
% imgParam      : structure array contains all image data
%
% Output
% --------------
% DIMWIParams   : algorithm parameters for DIMWI
% icvf          : intra-cellular volume fraction
% ff            : fibre fraction
% theta         : fibre orientation with respect to B0 direction
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 10 August 2021
% Date modified:
%
%
function [DIMWIParams, icvf, ff, theta] = setup_DIMWI(imgParam, algoParam)

DIMWIParams = algoParam.DIMWI;
b0dir       = imgParam.b0dir;

dims	= size(imgParam.img);   
dims    = dims(1:3);

% check if intra-axonal water volume fraction is needed
if DIMWIParams.isVic
    icvf = double(imgParam.icvf);
else
    icvf = zeros(dims);
end

% check if fibre orientation is needed
if DIMWIParams.isFreqMW || DIMWIParams.isFreqIW || DIMWIParams.isR2sEW
    ff    = double(imgParam.ff); % fibre fraction
    try 
        fo    = double(imgParam.fo); % fibre orientation w.r.t. B0
        theta = zeros(size(ff));
        for kfo = 1:size(fo,5)
            theta(:,:,:,kfo) = AngleBetweenV1MapAndB0(fo(:,:,:,:,kfo),b0dir);
        end
    catch 
        theta = double(imgParam.theta); % theta map
    end
    % normalise fibre fraction
    ff = bsxfun(@rdivide,ff,sum(ff,4));
    ff(isnan(ff)) = 0;
else
    ff      = ones(dims);
    theta   = zeros(dims);
end

% store fixed parameters
fixParam.b0     = double(imgParam.b0);
fixParam.rho_mw = double(imgParam.rho_mw);
fixParam.E      = double(imgParam.E);
fixParam.x_i    = double(imgParam.x_i);
fixParam.x_a    = double(imgParam.x_a);
DIMWIParams.fixParam = fixParam;

end