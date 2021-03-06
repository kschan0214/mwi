%% fitRes = mwi_3cx_T2s(algoPara,imgPara)
%
% Input
% --------------
% algoPara.maxIter : maximum iteration allows (default 500)
% algoPara.isROI   : boolean ROI analysis (default false)
% algoPara.DEBUG   : debug mode (default false)
% algoPara.fcnTol  : function tolerance (default: 1e-5)
% algoPara.isWeighted : boolean cost function weighted by echo intensity (default: True)
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
% Date last modified: 16 August 2018
%
%
function fitRes = mwi_3cx_T2s_HCFM_MWFreq_multiTol(algoPara,imgPara)
disp('Myelin water imaing: ME-T2* model');
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

% capture all algorithm parameters
numMagn    = algoPara.numMagn;
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
isWeighted = algoPara.isWeighted;
isParallel = algoPara.isParallel;
userDefine = algoPara.userDefine;
isInvivo   = algoPara.isInvivo;
% isROI      = algoPara.isROI; 

% if DEBUG is on then disables parallel computing
if DEBUG
    isParallel = false;
end

te    = imgPara.te;
data  = imgPara.img;
mask  = imgPara.mask;
fm    = imgPara.fieldmap;
pini  = imgPara.pini;

% for complex fitting we doubled the elements in the cost function
% fcnTol = fcnTol./(2*numel(te)-numMagn);
% fcnTol = fcnTol./(numel(te));
% if numMagn~=numel(te)
%     fcnTol = fcnTol*10;
% end

% display fitting message
if verbose
    disp('The following fitting parameters are used:');
    fprintf('Max. iterations = %i\n',maxIter);
    fprintf('Function tolerance = %e\n',fcnTol);
    if isWeighted
        disp('Cost function weighted by echo intensity: True');
    else
        disp('Cost function weighted by echo intensity: False');
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
        
    % initial guess for in-vivo case
    if isInvivo
        disp('Initial guesses for in vivo study');
    else
        disp('Initial guesses for ex vivo study');
    end
end

[ny,nx,nz,~,~] = size(data);

options_step1 = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*15,...
    'StepTolerance',1e-6,'FunctionTolerance',1e-5);

options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*15,...
    'StepTolerance',1e-6,'FunctionTolerance',fcnTol);

% get the order of magnitude of signal
% orderS = floor(log10(max(abs(data(:)))));
% fdss = [10^(orderS-6),10^(orderS-6),10^(orderS-6),1e-8,1e-8,1e-8,1e-8,1e-8,1e-8],1e-8]];
% options.FiniteDifferenceStepSize = fdss;
if ~DEBUG
    options.Display = 'off';
    options_step1.Display = 'off';
end

if numMagn==numel(te)   % magnitude fitting has 8 estimates
    estimates_step1 = zeros(ny,nx,nz,8);
else
    estimates_step1 = zeros(ny,nx,nz,8); % others have 10 estimates
end
if numMagn==numel(te)   % magnitude fitting has 8 estimates
    estimates_step2 = zeros(ny,nx,nz,6);
else
    estimates_step2 = zeros(ny,nx,nz,6); % others have 10 estimates
end

resnorm   = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        if verbose
            fprintf('Processing slice %i\n',kz);
        end
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % T2*w
                    s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0 = fm(ky,kx,kz);
                    pini0 = pini(ky,kx,kz);
                    [estimates_step1(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,pini0,numMagn,isWeighted,userDefine,isInvivo,options_step1,DEBUG);
                    [estimates_step2(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel_step2(s,te,db0,pini0,squeeze(estimates_step1(ky,kx,kz,:)),numMagn,isWeighted,userDefine,isInvivo,options,DEBUG);
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
                    % T2*w
                    s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0 = fm(ky,kx,kz);
                    pini0 = pini(ky,kx,kz);
                    [estimates_step1(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,pini0,numMagn,isWeighted,userDefine,isInvivo,options_step1,DEBUG);
                    [estimates_step2(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel_step2(s,te,db0,pini0,squeeze(estimates_step1(ky,kx,kz,:)),numMagn,isWeighted,userDefine,isInvivo,options,DEBUG);
                end
            end
        end
    end
end

fitRes.estimates_step1 = estimates_step1;
fitRes.estimates_step2 = estimates_step2;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,te,db0,pini0,numMagn,isWeighted,userDefine,isInvivo,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end

% [~,~,s0] = R2star_regression(s,te);
b = [ones(length(te),1), -te(:)]\log(abs(s(:)));
% r2s(kx,ky,kz) = b(2);
s0 = exp(b(1));
if s0<0
    s0=0;
end
Amw0 = 0.1*abs(s(1));   Amwlb = 0;  Amwub = 2*abs(s0);
Aiw0 = 0.6*abs(s(1));   Aiwlb = 0;  Aiwub = 2*abs(s0);
Aew0 = 0.3*abs(s(1));	Aewlb = 0;	Aewub = 2*abs(s0);        
% set initial guesses
if isInvivo
    % in vivo reference
    t2smy0 = 10;            t2smylb = 1;     t2smyub = 25;
    t2sax0 = 64;           	t2saxlb = 25;    t2saxub = 150;
    t2sex0 = 48;           	t2sexlb = 25;    t2sexub = 150;
else
    % ex vivo reference
    t2smy0 = 10;            t2smylb = 1;     t2smyub = 25;
    t2sax0 = 54;           	t2saxlb = 25;    t2saxub = 150;
    t2sex0 = 38;           	t2sexlb = 25;    t2sexub = 150;
end
    theta0 = pi/4;      thetalb = -pi; theta0ub = pi;

if numMagn==numel(te) % magnitude fitting
    fmy0   = 5;           	fmylb   = 5-75;  	fmyub   = 5+75;
    fbkg0   = 0;          	fbkglb   = -25;      fbkgub   = +25;

    % set initial guess and fitting boundaries
    x0 = double([Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,fmy0,fax0]);
    lb = double([Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,fmylb,faxlb]);
    ub = double([Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,fmyub,faxub]);
else    % other fittings
    fmy0   = 5;          	fmylb   = -75; 	fmyub   = +75;
    fbkg0   = db0;          	fbkglb   = db0-25;  	fbkgub   = db0+25;
    if isnan(pini0)
        pini0  = angle(exp(1i*(-2*pi*db0*te(1)-angle(s(1)))));
    end
    pinilb = -2*pi;         piniub=2*pi;

    % set initial guess and fitting boundaries
%     x0 = double([Amw0,Aiw0,Aew0,t2smy0,t2sax0,fmy0,fbkg0,pini0]);
%     lb = double([Amwlb,Aiwlb,Aewlb,t2smylb,t2saxlb,fmylb,fbkglb,pinilb]);
%     ub = double([Amwub,Aiwub,Aewub,t2smyub,t2saxub,fmyub,fbkgub,piniub]);
    x0 = double([Amw0,Aiw0,Aew0,t2smy0,t2sax0,theta0,fbkg0,pini0]);
    lb = double([Amwlb,Aiwlb,Aewlb,t2smylb,t2saxlb,thetalb,fbkglb,pinilb]);
    ub = double([Amwub,Aiwub,Aewub,t2smyub,t2saxub,theta0ub,fbkgub,piniub]);
end

% set initial guess and fitting bounds here
if ~isempty(userDefine.x0)
%     x0 = userDefine.x0;
    x0(~isnan(userDefine.x0)) = userDefine.x0(~isnan(userDefine.x0));
end
if ~isempty(userDefine.lb)
%     lb = userDefine.lb;
    lb(~isnan(userDefine.lb)) = userDefine.lb(~isnan(userDefine.lb));
end
if ~isempty(userDefine.ub)
%     ub = userDefine.ub;
    ub(~isnan(userDefine.ub)) = userDefine.ub(~isnan(userDefine.ub));
end

% run fitting algorithm here
[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,numMagn,isWeighted,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,numMagn,isWeighted,DEBUG)
% distribute fitting parameters
A_mw=x(1);    A_iw=x(2);   A_ew=x(3);
t2s_mw=x(4)*1e-3;  t2_fw=x(5)*1e-3; 
theta=x(6); 

if numMagn==numel(te) % magnitude fitting
    freq_bkg=0;        pini=0;
else    % other fittings
    freq_bkg=x(7);     pini=x(8);
end

param.b0 = 3;
param.pini = pini;
param.freq_bkg = freq_bkg;

sHat = mwi_model_ssSPGR_3T2scc_HCFM(te,A_mw,A_iw,A_ew,...
                                            t2s_mw,t2_fw,...
                                            theta,param);

if size(sHat,1) ~= size(s,1)
    sHat = sHat.';
end

% compute fitting residual
if isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
    w = sqrt(abs(s)/norm(abs(s(:))));
else
    % compute the cost without weights (=same weights)
    w = sqrt(ones(s)/numel(s));
end
err = computeFiter(s,sHat,numMagn,w);

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
err = err ./ norm(abs(s(:)));

%%%%%%% TEST 20180307 %%%%%%%%
% if numMagn==0
%     err = err/2;
% end

% Debug module
if DEBUG
    Debug_display(s,sHat,err,te,x,numMagn);
end

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel_step2(s,te,db0,pini0,estimates,numMagn,isWeighted,userDefine,isInvivo,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end

% [~,~,s0] = R2star_regression(s,te);
b = [ones(length(te),1), -te(:)]\log(abs(s(:)));
% r2s(kx,ky,kz) = b(2);
s0 = exp(b(1));
if s0<0
    s0=0;
end
Amw0 = estimates(1);   Amwlb = 0;  Amwub = 2*abs(s0);
Aiw0 = estimates(2);   Aiwlb = 0;  Aiwub = 2*abs(s0);
Aew0 = estimates(3);	Aewlb = 0;	Aewub = 2*abs(s0);  

t2 = estimates(4:5);

theta0 = estimates(6);      thetalb = -pi; theta0ub = pi;

if numMagn==numel(te) % magnitude fitting
    fmy0   = 5;           	fmylb   = 5-75;  	fmyub   = 5+75;
    fbkg0   = 0;          	fbkglb   = -25;      fbkgub   = +25;

    % set initial guess and fitting boundaries
    x0 = double([Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,fmy0,fax0]);
    lb = double([Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,fmylb,faxlb]);
    ub = double([Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,fmyub,faxub]);
else    % other fittings
    fbkg0   = estimates(7);          	fbkglb   = db0-25;  	fbkgub   = db0+25;
    pini0  = estimates(8); pinilb = -2*pi;         piniub=2*pi;

    % set initial guess and fitting boundaries
    x0 = double([Amw0,Aiw0,Aew0,theta0,fbkg0,pini0]);
    lb = double([Amwlb,Aiwlb,Aewlb,thetalb,fbkglb,pinilb]);
    ub = double([Amwub,Aiwub,Aewub,theta0ub,fbkgub,piniub]);
end

% set initial guess and fitting bounds here
if ~isempty(userDefine.x0)
%     x0 = userDefine.x0;
    x0(~isnan(userDefine.x0)) = userDefine.x0(~isnan(userDefine.x0));
end
if ~isempty(userDefine.lb)
%     lb = userDefine.lb;
    lb(~isnan(userDefine.lb)) = userDefine.lb(~isnan(userDefine.lb));
end
if ~isempty(userDefine.ub)
%     ub = userDefine.ub;
    ub(~isnan(userDefine.ub)) = userDefine.ub(~isnan(userDefine.ub));
end

% run fitting algorithm here
[x,res] = lsqnonlin(@(y)CostFunc_step2(y,s,te,t2,numMagn,isWeighted,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc_step2(x,s,te,t2,numMagn,isWeighted,DEBUG)
% distribute fitting parameters
A_mw=x(1);    A_iw=x(2);   A_ew=x(3);
t2s_mw=t2(1)*1e-3;  t2_fw=t2(2)*1e-3; 
theta=x(4); 

if numMagn==numel(te) % magnitude fitting
    freq_bkg=0;        pini=0;
else    % other fittings
    freq_bkg=x(5);     pini=x(6);
end

param.b0 = 3;
param.pini = pini;
param.freq_bkg = freq_bkg;

sHat = mwi_model_ssSPGR_3T2scc_HCFM(te,A_mw,A_iw,A_ew,...
                                            t2s_mw,t2_fw,...
                                            theta,param);

if size(sHat,1) ~= size(s,1)
    sHat = sHat.';
end

% compute fitting residual
if isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
    w = sqrt(abs(s)/norm(abs(s(:))));
else
    % compute the cost without weights (=same weights)
    w = sqrt(ones(s)/numel(s));
end
err = computeFiter(s,sHat,numMagn,w);

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
err = err ./ norm(abs(s(:)));

%%%%%%% TEST 20180307 %%%%%%%%
% if numMagn==0
%     err = err/2;
% end

% Debug module
if DEBUG
    Debug_display(s,sHat,err,te,x,numMagn);
end

end






%% check and set default
function [algoPara2,imgPara2,isValid]=CheckAndSetPara(algoPara,imgPara)

imgPara2 = imgPara;
algoPara2 = algoPara;
isValid = true;

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
    algoPara2.fcnTol = 1e-6;
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
% check # of phase-corrupted echoes
try
    algoPara2.isInvivo = algoPara.isInvivo;
catch
    algoPara2.isInvivo = true;
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
% check field map
try
    imgPara2.fieldmap = imgPara.fieldmap;
    disp('Field map input: True');
catch
    imgPara2.fieldmap = zeros(size(imgPara2.mask));
    disp('Field map input: False');
end
% check field map
try
    imgPara2.pini = imgPara.pini;
    disp('Initial phase input: True');
catch
    imgPara2.pini = ones(size(imgPara2.mask))*nan;
    disp('Initial phase input: False');
end

end

%% Info display for debug mode
function Debug_display(s,sHat,err,te,x,numMagn)
    global DEBUG_resnormAll
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    figure(99);

    if numMagn==numel(te)
        subplot(211);plot(te(:).',abs(permute(s(:),[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))*1.1]);
        title('Magnitude');
        plot(te(:).',abs(permute(sHat(:),[2 1])),'x-');plot(te(:).',(abs(permute(sHat(:),[2 1]))-abs(permute(s(:),[2 1]))),'ro-.');
        hold off;
        text(te(1)*0.5,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
%         text(te(1)*0.5,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f',...
%             Amy,Aax,Aex,t2smy,t2sax,t2sex,fmybg,faxbg));
        subplot(212);plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm')
    else
        subplot(411);
        plot(te(:).',real(permute(s,[2 1])),'k^-');hold on;
        plot(te(:).',imag(permute(s,[2 1])),'ks-');
        ylim([min([real(s(:));imag(s(:))])*1.1,max([real(s(:));imag(s(:))])*1.1]);
        title('Real and Imaginary');
        subplot(412);plot(te(:).',abs(permute(s,[2 1])),'k^-');hold on;
        ylim([0,max(abs(s(:)))*1.1]);
        title('Magnitude');
        subplot(413);plot(te(:).',angle(permute(s,[2 1])),'k^-');hold on;
        ylim([-4 4]);
        title('Phase');
        subplot(411);
        plot(te(:).',real(permute(sHat(:),[2 1])),'bx-');
        plot(te(:).',imag(permute(sHat(:),[2 1])),'b*-');
        plot(te(:).',real(permute(sHat(:),[2 1]))-real(permute(s(:),[2 1])),'rx-.');
        plot(te(:).',imag(permute(sHat(:),[2 1]))-imag(permute(s(:),[2 1])),'r*-.');hold off;
        subplot(412);
        plot(te(:).',abs(permute(sHat(:),[2 1])),'x-');
        plot(te(:).',abs(abs(permute(sHat(:),[2 1]))-abs(permute(s(:),[2 1]))),'ro-.');hold off;
%         text(te(1)*0.5,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
%         text(te(1)*0.5,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f,fex=%f,pini=%f',...
%             Amy,Aax,Aex,t2smy,t2sax,t2sex,fmybg,faxbg,fexbg,pini));
        subplot(413);
        plot(te(:).',angle(permute(sHat(:),[2 1])),'x-');
        plot(te(:).',angle(permute(s(:).*conj(sHat(:)),[2 1])),'ro-.');hold off;
        subplot(414);
        plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm');
    end
    if length(DEBUG_resnormAll) <300
        xlim([0 300]);
    else
        xlim([length(DEBUG_resnormAll)-300 length(DEBUG_resnormAll)]);
    end
    drawnow;
end