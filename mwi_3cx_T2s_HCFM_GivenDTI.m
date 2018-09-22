%% fitRes = mwi_3cx_T2s_wharton(algoPara,imgPara)
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
% Date created: 23 May 2018
% Date last modified: 
%
%
function fitRes = mwi_3cx_T2s_HCFM_GivenDTI(algoPara,imgPara)
disp('Myelin water imaing: ME-T2* Wharton and Bowtell model');
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
numMagn    = algoPara.numMagn;
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

te    = imgPara.te;
data  = imgPara.img;
mask  = imgPara.mask;
fm    = imgPara.fieldmap;
b0    = imgPara.b0;
theta = imgPara.theta;

% for complex fitting we doubled the elements in the cost function
% fcnTol = fcnTol./(2*numel(te)-numMagn);
% fcnTol = fcnTol./(numel(te));
if numMagn~=numel(te)
    fcnTol = fcnTol*10;
end

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

if numMagn==numel(te)   % magnitude fitting has 8 estimates
    estimates = zeros(ny,nx,nz,8);
else
    estimates = zeros(ny,nx,nz,10); % others have 10 estimates
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
                    theta0 = theta(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,theta0,b0,db0,numMagn,isWeighted,userDefine,options,DEBUG);
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
                    theta0 = theta(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,theta0,b0,db0,numMagn,isWeighted,userDefine,options,DEBUG);
                end
            end
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,te,theta,b0,db0,numMagn,isWeighted,userDefine,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end

b = [ones(length(te),1), -te(:)]\log(abs(s(:)));
% r2s(kx,ky,kz) = b(2);
s0 = exp(b(1));

fvf0    = 0.5;      fvflb   = 0;        fvfub   = 1;
g0      = 0.8;      glb     = 0;        gub     = 1;
pd_m0   = 0.7;    	pd_mlb  = 0;        pd_mub  = 1;
% e0      = 0.02;     elb     = 0.02;     eub     = 0.02;
t2_m0   = 15e-3;    t2_mlb  = 1e-3;     t2_mub  = 25e-3;
t2_n0   = 64e-3;    t2_nlb  = 25e-3;    t2_nub  = 150e-3;
chi_i0  = 0;     chi_ilb = -0.2;  	chi_iub = 0.2;
chi_a0  = 0;     chi_alb = -0.2;     chi_aub = 0.2;
% theta0  = 45;       thetalb = 0;        thetaub = 90; 
if s0>0
s00 = s0; s0lb = 0; s0ub=2*s0;
else
    s00 = 0; s0lb = 0; s0ub=2*s0;
end


% set initial guess and fitting boundaries
% x0 = double([fvf0, g0, pd_m0, e0, t2_m0, t2_n0, chi_i0, chi_a0, theta0]);
% lb = double([fvflb,glb,pd_mlb,elb,t2_mlb,t2_nlb,chi_ilb,chi_alb,thetalb]);
% ub = double([fvfub,gub,pd_mub,eub,t2_mub,t2_nub,chi_iub,chi_aub,thetaub]);
% x0 = double([fvf0, g0, pd_m0, t2_m0, t2_n0, chi_i0, chi_a0]);
% lb = double([fvflb,glb,pd_mlb,t2_mlb,t2_nlb,chi_ilb,chi_alb]);
% ub = double([fvfub,gub,pd_mub,t2_mub,t2_nub,chi_iub,chi_aub]);
x0 = double([fvf0, g0, pd_m0, t2_m0, t2_n0, chi_i0, chi_a0,s00]);
lb = double([fvflb,glb,pd_mlb,t2_mlb,t2_nlb,chi_ilb,chi_alb,s0lb]);
ub = double([fvfub,gub,pd_mub,t2_mub,t2_nub,chi_iub,chi_aub,s0ub]);

if numMagn~=numel(te) % magnitude fitting
 
    f_bkg0	= db0;          	f_bkglb   = db0-25;  	f_bkgub   = db0+25;
    pini0	= angle(exp(1i*(-2*pi*db0*te(1)-angle(s(1)))));        pinilb = -pi;         piniub=pi;
    
    % set initial guess and fitting boundaries
    x0 = double([x0,f_bkg0,pini0]);
    lb = double([lb,f_bkglb,pinilb]);
    ub = double([ub,f_bkgub,piniub]);
    
end

% set initial guess and fitting bounds if they are provided
if ~isempty(userDefine.x0)
    x0(~isnan(userDefine.x0)) = userDefine.x0(~isnan(userDefine.x0));
end
if ~isempty(userDefine.lb)
    lb(~isnan(userDefine.lb)) = userDefine.lb(~isnan(userDefine.lb));
end
if ~isempty(userDefine.ub)
    ub(~isnan(userDefine.ub)) = userDefine.ub(~isnan(userDefine.ub));
end

% run fitting algorithm here
[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,theta,numMagn,isWeighted,b0,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,theta,numMagn,isWeighted,b0,DEBUG)
% obatin fitting parameters
fvf     = x(1);    
g       = x(2);   
pd_m	= x(3);
% e       = x(4);
e       = 0.02;  
t2_m    = x(4); t2_n    = x(5);
chi_i   = x(6); chi_a   = x(7);
s0 = x(8);
% theta   = x(9);

if numMagn==numel(te) % magnitude fitting
    f_bkg=0;        pini=0;
else    % other fittings
%     f_bkg=x(8);    pini=x(9);
f_bkg=x(9);    pini=x(10);
end

% simulate signal based on parameter input
sHat = mwi_model_3cc_HCFM(te, ...
                          fvf, g, pd_m, e, ...
                          t2_m, t2_n, chi_i, chi_a, theta,...
                          f_bkg, pini, s0, b0);

% s       = s     ./ norm(abs(s));
% sHat    = sHat  ./ norm(abs(sHat));

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

%%%%%%% TEST 20180307 %%%%%%%%
% if numMagn==0
%     err = err/2;
% end

% Debug module
if DEBUG
    global DEBUG_resnormAll
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    figure(99);
    if numMagn==numel(te)
        subplot(211);plot(te(:).',abs(permute(s(:),[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))*1.1]);
        title('Magnitude');
        plot(te(:).',abs(permute(sHat(:),[2 1])),'x-');plot(te(:).',(abs(permute(sHat(:),[2 1]))-abs(permute(s(:),[2 1]))),'ro-.');
        hold off;
        text(te(1)*0.5,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
        text(te(1)*0.5,max(abs(s(:))*0.1),sprintf('fvf=%f,g=%f,pd_m=%f,e=%f,t2_m=%f,t2_n=%f,fmy=%f,fax=%f,theta=%f',...
            fvf,g,pd_m,e,t2_m,t2_n,chi_i,chi_a,theta));
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
        plot(te(:).',real(permute(sHat,[2 1])),'bx-');
        plot(te(:).',imag(permute(sHat,[2 1])),'b*-');
        plot(te(:).',real(permute(sHat,[2 1]))-real(permute(s,[2 1])),'rx-.');
        plot(te(:).',imag(permute(sHat,[2 1]))-imag(permute(s,[2 1])),'r*-.');hold off;
        subplot(412);
        plot(te(:).',abs(permute(sHat,[2 1])),'x-');
        plot(te(:).',abs(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'ro-.');hold off;
        text(te(1)*0.5,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
        text(te(1)*0.5,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f,fex=%f,pini=%f',...
            Amy,Aax,Aex,t2smy,t2sax,t2sex,fmybg,faxbg,fexbg,pini));
        subplot(413);
        plot(te(:).',angle(permute(sHat,[2 1])),'x-');
        plot(te(:).',angle(permute(s.*conj(sHat),[2 1])),'ro-.');hold off;
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
% check B0
try
    imgPara2.b0 = imgPara.b0;
    disp('Field strength input: True');
catch
    imgPara2.b0 = 3;
    disp('Field strength input: False');
end

end