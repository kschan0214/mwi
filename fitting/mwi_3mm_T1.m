%% fitRes = mwi_3cx_T2s(algoPara,imgPara)
%
% Input
% --------------
% algoPara.maxIter : maximum iteration allows (default 500)
% algoPara.isROI   : boolean ROI analysis (default false)
% algoPara.DEBUG   : debug mode (default false)
% algoPara.fcnTol  : function tolerance (default: 1e-5)
% imgPara.img      : 4D image data, time in 4th dimension
% imgPara.mask     : signal mask
% imgPara.te       : echo times
% imgPara.fieldmap : background field (default: 0)
%
% Output
% --------------
% fitRes.estimates : fitting estimates (Ampl_n,t2s_n,freq_n)
% fitres.resnorm   : L2 norm of fitting residual
%
% Description: Myelin water mapping by fitting complex mode(c) with
% complex-valued data(c) or magnitude data(m) (ref. Model 3 in Nam et al. 2015 NeuroImage)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 February 2018
% Date last modified: 
%
%
function fitRes = mwi_3mm_T1(algoPara,imgPara)
% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

% get debug mode and verbose
DEBUG   = algoPara.DEBUG;
verbose = algoPara.verbose;

% capture all parameters
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
isWeighted = algoPara.isWeighted;
% isROI      = algoPara.isROI; 
isParallel = algoPara.isParallel;
userDefine = algoPara.userDefine;

% if DEBUG on disables parallel computing
if DEBUG
    isParallel = false;
end

fa    = imgPara.fa;
tr    = imgPara.tr;
data  = abs(imgPara.img);
mask  = imgPara.mask;
b1map = imgPara.b1map;

numMagn    = length(fa);

% for complex fitting we doubled the elements in the cost function
% fcnTol = fcnTol./(2*numel(te)-numMagn);
fcnTol = fcnTol./(numel(fa));

% display fitting message
if verbose
    disp('The following fitting parameters are used:');
    fprintf('Max. iterations = %i\n',maxIter);
    fprintf('Function tolerance = %e\n',fcnTol);
    if isWeighted
        disp('Cost function weighted by signal intensity: True');
    else
        disp('Cost function weighted by signal intensity: False');
    end
    % initial guess and fitting bounds
    if isempty(userDefine.x0)
        disp('Default initial guess: True');
    else
        disp('Default initial guess: False');
    end
    if isempty(userDefine.lb)
        disp('Default lower bound: True');
    else
        disp('Default lower bound: False');
    end
    if isempty(userDefine.ub)
        disp('Default upper bound: True');
    else
        disp('Default upper bound: False');
    end
end

[ny,nx,nz,~,~] = size(data);

options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*15,...
    'StepTolerance',1e-6,'FunctionTolerance',fcnTol);

% get the order of magnitude of signal
% orderS = floor(log10(max(abs(data(:)))));
% fdss = [10^(orderS-6),10^(orderS-6),10^(orderS-6),1e-8,1e-8,1e-8,1e-8,1e-8,1e-8],1e-8]];
% options.FiniteDifferenceStepSize = fdss;
if ~DEBUG
    options.Display = 'off';
end

% if numMagn==numel(te)   % magnitude fitting has 8 estimates
    estimates = zeros(ny,nx,nz,5);
% end

resnorm   = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % T1w
                    s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    b1 = b1map(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,tr,b1,numMagn,isWeighted,userDefine,options,DEBUG);
                end
            end
        end
    end
else
    for kz=1:nz
        for ky=1:ny
            for kx=1:nx
                if mask(ky,kx,kz)>0
                    % T1w
                    s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    b1 = b1map(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,tr,b1,numMagn,isWeighted,userDefine,options,DEBUG);
                end
            end
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,tr,b1,numMagn,isWeighted,userDefine,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end
% SMAX = max(s(:));
% [~,m0] = DESPOT1(s,fa,tr);
B_t1 = mwi_GetSpeciesT1Decay(fa,tr,[537e-3,1183e-3,1183e-3]);
A_t1 = pinv(B_t1)*abs(s(:));

% 
% Amy0   = 0.1*abs(m0);   	Amylb   = 0;        Amyub   = 1*abs(m0);
% Aax0   = 0.6*abs(m0);       Aaxlb   = 0;        Aaxub   = 1*abs(m0);
% Aex0   = 0.3*abs(m0);    	Aexlb   = 0;        Aexub   = 1*abs(m0);
t1my0  = 537e-3;           	t1mylb  = 100e-3;   t1myub  = 650e-3;
t1ax0  = 1183e-3;          	t1axlb  = 650e-3;   t1axub  = 2000e-3;
% t1ex0 = 48e-3;           	t1exlb = 25e-3;     t1exub = 150e-3;

Amy0   = A_t1(1);   	Amylb   = 0;        Amyub   = 1*abs(Amy0);
Aax0   = A_t1(2);       Aaxlb   = 0;        Aaxub   = 2*abs(Aax0);
Aex0   = A_t1(3);    	Aexlb   = 0;        Aexub   = 2*abs(Aex0);


    x0 = double([Amy0 ,Aax0 ,Aex0 ,t1my0 ,t1ax0]);
    lb = double([Amylb,Aaxlb,Aexlb,t1mylb,t1axlb]);
    ub = double([Amyub,Aaxub,Aexub,t1myub,t1axub]);


% set initial guess and fitting bounds if they are provided
if ~isempty(userDefine.x0)
    x0 = userDefine.x0;
end
if ~isempty(userDefine.lb)
    lb = userDefine.lb;
end
if ~isempty(userDefine.ub)
    ub = userDefine.ub;
end

% run fitting algorithm here
[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,tr,b1,numMagn,isWeighted,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,tr,b1,numMagn,isWeighted,DEBUG)
% obatin fitting parameters
Amy=x(1);    Aax=x(2);   Aex=x(3);
t1my=x(4);  t1ax=x(5); 

% simulate signal based on parameter input
% if numMagn==numel(te)
sHat = mwi_model_3T1(fa,tr,Amy,Aax,Aex,t1my,t1ax,t1ax,b1);
% else
%     sHat = mwi_model_3cc_nam2015(te,Amy,Aax,Aex,t1my,t1ax,t2sex,fmybg,faxbg,fexbg,pini);
% end

% compute fitting residual
if isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
    w = abs(s(:))/norm(abs(s(:)));
    err = computeFiter(s,sHat,numMagn,w);
%     w = abs(s(:));
%     err = err.*w(:);
else
    err = computeFiter(s,sHat,numMagn);
end

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
err = err ./ norm(abs(s(:)));

% Debug module
if DEBUG
    global DEBUG_resnormAll
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    figure(99);
        subplot(211);plot(fa(:).',abs(permute(s(:),[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))*1.1]);
        title('Magnitude');
        plot(fa(:).',abs(permute(sHat(:),[2 1])),'x-');plot(fa(:).',(abs(permute(sHat(:),[2 1]))-abs(permute(s(:),[2 1]))),'ro-.');
        hold off;
        text(fa(1)*0.5,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
        text(fa(1)*0.5,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t1my=%f,t1ax=%f',...
            Amy,Aax,Aex,t1my,t1ax));
        subplot(212);plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm')
    if length(DEBUG_resnormAll) <300
        xlim([0 300]);
    else
        xlim([length(DEBUG_resnormAll)-300 length(DEBUG_resnormAll)]);
    end
    drawnow;
end

end

%% check and set default
function [algoPara2,imgPara2,isValid]=CheckAndSetPara(algoPara,imgPara)

imgPara2 = imgPara;
algoPara2 = algoPara;
isValid = true;

% check if the number of echo times matches with the data
if length(imgPara.fa) ~= size(imgPara.img,4)
    isValid = false;
end

% check debug
try
    algoPara2.DEBUG = algoPara.DEBUG;
catch
    algoPara2.DEBUG = false;
end
% check verbose
try
    algoPara2.verbose = algoPara.verbose;
catch
    algoPara2.verbose = true;
end

% check maximum iterations allowed
try
    algoPara2.maxIter = algoPara.maxIter;
catch
    algoPara2.maxIter = 500;
end
% check parallel computing 
try
    algoPara2.isParallel = algoPara.isParallel;
catch
    algoPara2.isParallel = false;
end
% check function tolerance
try
    algoPara2.fcnTol = algoPara.fcnTol;
catch
    algoPara2.fcnTol = 1e-5;
end
% check weighted sum of cost function
try
    algoPara2.isWeighted = algoPara.isWeighted;
catch
    algoPara2.isWeighted = true;
end

% check ROI analysis
try
    algoPara2.isROI = algoPara.isROI;
catch
    algoPara2.isROI = false;
end
% check # of phase-corrupted echoes
try
    algoPara2.numMagn = algoPara.numMagn;
catch
    algoPara2.numMagn = numel(imgPara.te);
end
% check user bounds and initial guesses
try
    algoPara2.userDefine.x0 = algoPara.userDefine.x0;
catch
    algoPara2.userDefine.x0 = [];
end
try
    algoPara2.userDefine.lb = algoPara.userDefine.lb;
catch
    algoPara2.userDefine.lb = [];
end
try
    algoPara2.userDefine.ub = algoPara.userDefine.ub;
catch
    algoPara2.userDefine.ub = [];
end

% check signal mask
try
    imgPara2.mask = imgPara.mask;
    disp('Default mask: False');
catch
    imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
    disp('Default mask: True');
end

% check b1 map
try
    imgPara2.b1map = imgPara.b1map;
    disp('Default B1: False');
catch
    imgPara2.b1map = ones(size(imgPara.img,1),size(imgPara.img,2),size(imgPara.img,3));
    disp('Default B1: True');
end

end