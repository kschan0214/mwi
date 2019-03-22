%% fitRes = mwi_2T13T2s_cm_VFAT2s(algoPara,imgPara)
%
% Input
% --------------
% algoPara.maxIter : maximum iteration allows (default 500)
% algoPara.DEBUG   : debug mode (default false)
% algoPara.fcnTol  : function tolerance (default: 1e-5)
% algoPara.isWeighted : boolean cost function weighted by echo intensity (default: True)
% imgPara.img      : 5D image data, time in 4th dimension, T1w in 5th dim
% imgPara.mask     : signal mask
% imgPara.te       : echo times
% imgPara.fa       : flip angles
% imgPara.b1map    : B1 map
%
% Output
% --------------
% fitRes.estimates : fitting estimates (Ampl_n,t2s_n,t1_n,freq_n)
% fitres.resnorm   : L2 norm of fitting residual
%
% Description: Myelin water mapping by fitting complex mode(c) with
% complex-valued data(c) or magnitude data(m) (ref. Model 3 in Nam et al. 2015 NeuroImage)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 January 2018
% Date last modified: 08 march 2018 
%
%
function fitRes = mwi_2T13T2s_cm_VFAT2s(algoPara,imgPara)
disp('Myelin water imaing: VFA-ME-T2* model');
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

% capture all fitting settings
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
% numMagn    = algoPara.numMagn;
isWeighted = algoPara.isWeighted;
weightMethod = algoPara.weightMethod;
userDefine = algoPara.userDefine;
npulse = algoPara.npulse;
isParallel = algoPara.isParallel;

if strcmpi(weightMethod,'norm') && strcmpi(weightMethod,'1stEcho')
    disp('Do not support the input weighting method');
    weightMethod = 'norm';
end

% % if DEBUG on disables parallel computing
% if DEBUG
%     isParallel = false;
% end

% capture all images related data
te    = imgPara.te;
tr    = imgPara.tr;
fa    = imgPara.fa;
data  = imgPara.img;
mask  = imgPara.mask;
b1map = imgPara.b1map;

% Magnitude fitting
numMagn    = length(imgPara.te);

% display fitting message
if verbose
    disp('The following fitting parameters are used:');
    fprintf('Max. iterations = %i\n',maxIter);
    fprintf('Function tolerance = %e\n',fcnTol);
    fprintf('No. of pulses for EPG-X = %i\n',npulse);
    if isWeighted
        disp('Weighted Cost function: True');
        disp(['Weighting method: ' weightMethod]);
    else
        disp('Weighted Cost function: False');
    end
    % type of fitting
    if numMagn==0
        disp('Fitting complex model with complex data');
    elseif numMagn==numel(te)
        disp('Fitting complex model with magnitude data');
    else
        fprintf('Fitting complex model with %i magnitude data and %i complex data\n',numMagn,numel(te)-numMagn);
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

% lsqnonlin setting
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*11,...
        'FunctionTolerance',fcnTol,'StepTolerance',1e-6);
  
% fdss = [1e-4,1e-4,1e-4,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8];  
% options.FiniteDifferenceStepSize = fdss;

% if DEBUG then display fitting message
if ~DEBUG
    options.Display = 'off';
end

% we are fitting 11 parameters with this model
estimates = zeros(ny,nx,nz,11);
resnorm   = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        if verbose
            fprintf('Processing slice %i\n',kz);
        end
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1 = b1map(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,npulse,numMagn,isWeighted,weightMethod,userDefine,options,DEBUG);
                end
            end
        end
    end
else
    for kz=1:nz
        if verbose
            fprintf('Processing slice %i\n',kz);
        end
        for ky=1:ny
            for kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1 = b1map(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,npulse,numMagn,isWeighted,weightMethod,userDefine,options,DEBUG);
                end
            end
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b1,npulse,numMagn,isWeighted,weightMethod,userDefine,options,DEBUG)
% define initial guesses
% estimate rho0 of the first echo
[~,rho0] = DESPOT1(abs(s(:,1)),fa,tr,'b1',b1);
% estimate t1 from later echo
[t10,~] = DESPOT1(abs(s(:,end-3)),fa,tr,'b1',b1);

% in case rho0 and t10 estimation go wrong then use defualt values
if rho0<0
    rho0=max(abs(s));
end
if t10<0
    t10=1000e-3;
end

Amy0   = 0.1*rho0;            Amylb   = 0;        Amyub   = 1*rho0;
Aax0   = 0.6*rho0;            Aaxlb   = 0;        Aaxub   = 1*rho0;
Aex0   = 0.3*rho0;            Aexlb   = 0;        Aexub   = 1*rho0;
t2smy0 = 10e-3;     t2smylb = 1e-3;     t2smyub = 25e-3;
t2sax0 = 64e-3; 	t2saxlb = 25e-3;    t2saxub = 200e-3;
t2sex0 = 48e-3; 	t2sexlb = 25e-3;    t2sexub = 200e-3;
fmy0   = 5;                fmylb   = 5-75;    fmyub   = 5+75;
fax0   = 0;                 faxlb   = -25;      faxub   = +25;
t1my0  = 300e-3;  	t1mylb  = 50e-3;    t1myub  = 650e-3;
t1l0   = t10;     	t1llb  = 500e-3;    t1lub  = 2000e-3;
kx0 = 2;    kxlb=0;   kxub = 6;

% set initial guess and fitting bounds here
x0 = [Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,t1my0,t1l0,fmy0,fax0,kx0];
lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,t1mylb,t1llb,fmylb,faxlb,kxlb];
ub = [Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,t1myub,t1lub,fmyub,faxub,kxub];

if ~isempty(userDefine.x0)
    x0(~isnan(userDefine.x0)) = userDefine.x0(~isnan(userDefine.x0));
end
if ~isempty(userDefine.lb)
    lb(~isnan(userDefine.lb)) = userDefine.lb(~isnan(userDefine.lb));
end
if ~isempty(userDefine.ub)
    ub(~isnan(userDefine.ub)) = userDefine.ub(~isnan(userDefine.ub));
end

if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
    x0
end

phiCycle = RF_phase_cycle(npulse,50);
for kfa=1:length(fa)
T3D_all{kfa} = PrecomputeT(phiCycle,d2r(fa(kfa)*b1));
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,npulse,T3D_all,numMagn,isWeighted,weightMethod,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,npulse,T3D_all,numMagn,isWeighted,weightMethod,DEBUG)
% capture all fitting parameters
Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5); t2sex=x(6);
t1my=x(7);  t1l=x(8);
fmy=x(9);  fax=x(10);
kx=x(11);


% sHat = mwi_model_2cm_jointT1T2s_epgX(fa,te,tr,A0,mwf,t2ss,t2sl,t1s,t1l,fs,0,b1,kx);
sHat = mwi_model_2T13T2scm_jointT1T2s_epgX(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1l,fmy,fax,b1,kx,npulse,T3D_all);

% err = computeFiter(s,sHat,numMagn);
% compute fitting residual
if isWeighted
    switch weightMethod
        case 'norm'
            % weighted the cost function by echo intensity, as suggested in Nam's paper
            w = abs(s(:))/norm(abs(s(:)));
        case '1stEcho'
            w = bsxfun(@rdivide,abs(s),abs(s(:,1)));
    end
    err = computeFiter(s,sHat,numMagn,w);
%     w = abs(s(:));
%     err = err.*w(:);
else
    err = computeFiter(s,sHat,numMagn);
end

% error weighted by echo strength
% w = abs(s(:))./norm(abs(s(:)));
% err = w.*err;

err = err./norm(abs(s));

% if DEBUG then plots current fitting result
if DEBUG
    global DEBUG_resnormAll
    figure(99);subplot(211);plot(te(:).',abs(permute(s,[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))+10]);
    title('Magnitude');
    plot(te(:).',abs(permute(sHat,[2 1])),'x-');plot(te(:).',(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'ro-.');
    hold off;
    text(te(1)/3,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
    text(te(1)/3,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f,T1my=%f,T1l=%f,kmy=%f',...
        Amy,Aax,Aex,t2smy,t2sax,t2sex,fmy,fax,t1my,t1l,kx));
    for kfa = 1:length(fa)
        text(te(1)/3,abs(s(kfa,1)),['FA ' num2str(fa(kfa))]);
    end
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    subplot(212);plot(DEBUG_resnormAll);
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

% check if the number of flip angles matches with the data
if length(imgPara.fa) ~= size(imgPara.img,5)
    isValid = false;
end
% check if the number of echo times matches with the data
if length(imgPara.te) ~= size(imgPara.img,4)
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
% check function tolerance
try
    algoPara2.fcnTol = algoPara.fcnTol;
catch
    algoPara2.fcnTol = 1e-5;
end
% check parallel computing 
try
    algoPara2.isParallel = algoPara.isParallel;
catch
    algoPara2.isParallel = false;
end
% check weighted sum of cost function
try
    algoPara2.isWeighted = algoPara.isWeighted;
catch
    algoPara2.isWeighted = true;
end
% check weighted sum of cost function
try
    algoPara2.weightMethod = algoPara.weightMethod;
catch
    algoPara2.weightMethod = 'norm';
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
    disp('Mask input: True');
catch
    imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
    disp('Mask input: false');
end

% check b1 map
try
    imgPara2.b1map = imgPara.b1map;
    disp('B1 input: True');
catch
    imgPara2.b1map = ones(size(imgPara.img,1),size(imgPara.img,2),size(imgPara.img,3));
    disp('B1 input: False');
end

end