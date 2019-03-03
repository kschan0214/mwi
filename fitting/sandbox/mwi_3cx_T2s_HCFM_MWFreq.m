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
function fitRes = mwi_3cx_T2s_HCFM_MWFreq(algoPara,imgPara)
disp('Myelin water imaing: ME-T2* model');
% check validity of the algorithm parameters and image parameters
[algoPara,imgPara]=CheckAndSetDefault(algoPara,imgPara);

% get debug mode 
DEBUG   = algoPara.DEBUG;

% capture all algorithm parameters
numMagn    = algoPara.numMagn;
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
isWeighted = algoPara.isWeighted;
isParallel = algoPara.isParallel;
userDefine = algoPara.userDefine;
isInvivo   = algoPara.isInvivo;

% if DEBUG is on then disables parallel computing
if DEBUG
    isParallel = false;
end

te    = imgPara.te;
data  = imgPara.img;
mask  = imgPara.mask;
fm    = imgPara.fieldmap;
pini  = imgPara.pini;

[ny,nx,nz,~,~] = size(data);

% set fitting options
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*15,...
    'StepTolerance',1e-6,'FunctionTolerance',fcnTol);

if ~DEBUG
    options.Display = 'off';
end

% create empty array for fitting results
if numMagn==numel(te)   % magnitude fitting has 6 estimates
    estimates = zeros(ny,nx,nz,6);
else
    estimates = zeros(ny,nx,nz,8); % others have 8 estimates
end

numMaskedVoxel = length(mask(mask==1));
numFittedVoxel = 0;
fprintf('%i voxles need to be fitted...\n',numMaskedVoxel);
progress='';

resnorm   = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % T2*w
                    s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0 = fm(ky,kx,kz);
                    pini0 = pini(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,pini0,numMagn,isWeighted,userDefine,isInvivo,options,DEBUG);
                end
            end
            % display progress
            tmp = mask(ky,:,kz);
            numFittedVoxel = numFittedVoxel + length(tmp(tmp==1));
            for ii=1:length(progress); fprintf('\b'); end
            progress=sprintf('Progress (%%): %d', floor(numFittedVoxel*100/numMaskedVoxel));
            fprintf(progress);
        end
    end
else
    for kz=1:nz
        fprintf('Processing slice %i\n',kz);
        for ky=1:ny
            for kx=1:nx
                if mask(ky,kx,kz)>0
                    % T2*w
                    s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0 = fm(ky,kx,kz);
                    pini0 = pini(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,pini0,numMagn,isWeighted,userDefine,isInvivo,options,DEBUG);
                end
            end
            % display progress
            tmp = mask(ky,:,kz);
            numFittedVoxel = numFittedVoxel + length(tmp(tmp==1));
            for ii=1:length(progress); fprintf('\b'); end
            progress=sprintf('Progress (%%): %d', floor(numFittedVoxel*100/numMaskedVoxel));
            fprintf(progress);
        end
    end
end
fprintf('\n');

fitRes.estimates = estimates;
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
      
% set initial guesses
if isInvivo
    % in vivo reference
    Amw0    = 0.1*abs(s(1));    Amwlb      = 0;    Amwub    = 2*abs(s0);
    Aiw0    = 0.6*abs(s(1));    Aiwlb      = 0;    Aiwub    = 2*abs(s0);
    Aew0    = 0.3*abs(s(1));    Aewlb      = 0;    Aewub    = 2*abs(s0);  
    t2smy0  = 10;               t2smylb    = 1;    t2smyub  = 25;
    t2fw0   = 64;               t2fwlb     = 25;   t2fwub   = 150;
else
    % ex vivo reference
    Amw0 = 0.2*abs(s(1));   Amwlb = 0;  Amwub = 2*abs(s0);
    Aiw0 = 0.6*abs(s(1));   Aiwlb = 0;  Aiwub = 2*abs(s0);
    Aew0 = 0.2*abs(s(1));	Aewlb = 0;	Aewub = 2*abs(s0);  
    t2smy0 = 10;            t2smylb = 1;     t2smyub = 25;
    t2fw0 = 54;           	t2fwlb = 25;    t2fwub = 150;
end
sin2theta0 = 0.5;      sin2thetalb = 0; sin2theta0ub = 1;

if numMagn==numel(te) % magnitude fitting
    % set initial guess and fitting boundaries
    x0 = double([Amy0,Aax0,Aex0,t2smy0,t2fw0,sin2theta0]);
    lb = double([Amylb,Aaxlb,Aexlb,t2smylb,t2fwlb,sin2thetalb]);
    ub = double([Amyub,Aaxub,Aexub,t2smyub,t2fwub,sin2theta0ub]);
    
else    % other fittings
    fbkg0   = db0;          	fbkglb   = db0-25;  	fbkgub   = db0+25;
    if isnan(pini0)
        pini0  = angle(exp(1i*(-2*pi*db0*te(1)-angle(s(1)))));
    end
    pinilb = -2*pi;         piniub=2*pi;

    % set initial guess and fitting boundaries
%     x0 = double([Amw0,Aiw0,Aew0,t2smy0,t2sax0,fmy0,fbkg0,pini0]);
%     lb = double([Amwlb,Aiwlb,Aewlb,t2smylb,t2saxlb,fmylb,fbkglb,pinilb]);
%     ub = double([Amwub,Aiwub,Aewub,t2smyub,t2saxub,fmyub,fbkgub,piniub]);
%     x0 = double([Amw0,Aiw0,Aew0,t2smy0,t2sax0,theta0,fbkg0,pini0]);
%     lb = double([Amwlb,Aiwlb,Aewlb,t2smylb,t2saxlb,thetalb,fbkglb,pinilb]);
%     ub = double([Amwub,Aiwub,Aewub,t2smyub,t2saxub,theta0ub,fbkgub,piniub]);
    x0 = double([Amw0,Aiw0,Aew0,t2smy0,t2fw0,sin2theta0,fbkg0,pini0]);
    lb = double([Amwlb,Aiwlb,Aewlb,t2smylb,t2fwlb,sin2thetalb,fbkglb,pinilb]);
    ub = double([Amwub,Aiwub,Aewub,t2smyub,t2fwub,sin2theta0ub,fbkgub,piniub]);
end

% set initial guess and fitting bounds here
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
[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,numMagn,isWeighted,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,numMagn,isWeighted,DEBUG)
% distribute fitting parameters
A_mw=x(1);    A_iw=x(2);   A_ew=x(3);
t2s_mw=x(4)*1e-3;  t2_fw=x(5)*1e-3; 
sin2theta=x(6); 

if numMagn==numel(te) 
    % magnitude fitting
    freq_bkg=0;        pini=0;
else
    % other fittings
    freq_bkg=x(7);     pini=x(8);
end

param.b0        = 3;
param.pini      = pini;
param.freq_bkg  = freq_bkg;

% simulate signal
sHat = mwi_model_ssSPGR_3T2scc_HCFM(te,A_mw,A_iw,A_ew,t2s_mw,t2_fw,sin2theta,param);

% matching the size of simulated signal and the measured signal
if size(sHat,1) ~= size(s,1)
    sHat = sHat.';
end

% compute fitting residual
if isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
    w = sqrt(abs(s)/norm(abs(s(:))));
else
    % compute the cost without weights (=same weights)
    w = sqrt(ones(size(s))/numel(s));
end
err = computeFiter(s,sHat,numMagn,w);

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
err = err ./ norm(abs(s(:)));

% Debug module
if DEBUG
    Debug_display(s,sHat,err,te,x,numMagn);
end

end

%% check and set default
function [algoPara2,imgPara2]=CheckAndSetDefault(algoPara,imgPara)

imgPara2 = imgPara;
algoPara2 = algoPara;

%%%%%%%%%% 1. check algorithm parameters %%%%%%%%%%
% check debug
try algoPara2.DEBUG = algoPara.DEBUG;                   catch; algoPara2.DEBUG = false; end
% check parallel computing 
try algoPara2.isParallel = algoPara.isParallel;         catch; algoPara2.isParallel = false; end
% check maximum iterations allowed
try algoPara2.maxIter = algoPara.maxIter;               catch; algoPara2.maxIter = 500; end
% check function tolerance
try algoPara2.fcnTol = algoPara.fcnTol;                 catch; algoPara2.fcnTol = 1e-6; end
% check weighted sum of cost function
try algoPara2.isWeighted = algoPara.isWeighted;         catch; algoPara2.isWeighted = true; end
% check # of phase-corrupted echoes
try algoPara2.numMagn = algoPara.numMagn;               catch; algoPara2.numMagn = numel(imgPara.te); end
% check # of phase-corrupted echoes
try algoPara2.isInvivo = algoPara.isInvivo;             catch; algoPara2.isInvivo = true; end
% check user bounds and initial guesses
try algoPara2.userDefine.x0 = algoPara.userDefine.x0;   catch; algoPara2.userDefine.x0 = [];end
try algoPara2.userDefine.lb = algoPara.userDefine.lb;   catch; algoPara2.userDefine.lb = [];end
try algoPara2.userDefine.ub = algoPara.userDefine.ub;   catch; algoPara2.userDefine.ub = [];end

%%%%%%%%%% 2. check data integrity %%%%%%%%%%
% check if the number of echo times matches with the data
if length(imgPara.te) ~= size(imgPara.img,4)
    error('The length of TE does not match with the last dimension of the image.');
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

%%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
disp('The following fitting parameters are used:');
fprintf('Max. iterations = %i\n',algoPara2.maxIter);
fprintf('Function tolerance = %e\n',algoPara2.fcnTol);
if algoPara2.isWeighted
    disp('Cost function weighted by echo intensity: True');
else
    disp('Cost function weighted by echo intensity: False');
end
% type of fitting
if algoPara2.numMagn==0
    disp('Fitting complex model with complex data');
elseif algoPara2.numMagn==numel(imgPara2.te)
    disp('Fitting complex model with magnitude data');
else
    fprintf('Fitting complex model with %i magnitude data and %i complex data\n',algoPara2.numMagn,numel(imgPara2.te)-algoPara2.numMagn);
end
% initial guess and fitting bounds
if isempty(algoPara2.userDefine.x0)
    disp('Default initial guess: True');
else
    disp('Default initial guess: False');
end
if isempty(algoPara2.userDefine.lb)
    disp('Default lower bound: True');
else
    disp('Default lower bound: False');
end
if isempty(algoPara2.userDefine.ub)
    disp('Default upper bound: True');
else
    disp('Default upper bound: False');
end
% initial guess for in-vivo case
if algoPara2.isInvivo
    disp('Initial guesses for in vivo study');
else
    disp('Initial guesses for ex vivo study');
end
    
end

%% Info display for debug mode
function Debug_display(s,sHat,err,te,x,numMagn)
    global DEBUG_resnormAll
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    figure(99);

    if numMagn==numel(te)
        subplot(2,2,1);plot(te(:).',abs(permute(s(:),[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))*1.1]); title('Magnitude');
        plot(te(:).',abs(permute(sHat(:),[2 1])),'x-.');hold off;
        subplot(2,2,2); plot(te(:).',(abs(permute(sHat(:),[2 1]))-abs(permute(s(:),[2 1]))),'ro-.'); title('residual');
        ha = subplot(2,2,3); pos = get(ha,'Position'); un = get(ha,'Units'); delete(ha)
        uitable('Data',x(:),'Units',un,'Position',pos);
        subplot(2,2,4);plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm')
        text(0.5,0.5,sprintf('resnorm=%f',sum(err(:).^2)),'Units','normalized');
    else
        subplot(2,3,1);
        plot(te(:).',abs(permute(s,[2 1])),'k^-');hold on;
        plot(te(:).',abs(permute(sHat(:),[2 1])),'x-.');hold off;
        ylim([min(abs(s(:)))*0.9,max(abs(s(:)))*1.1]);
        title('Magn.');
        
        subplot(2,3,2);
        plot(te(:).',angle(permute(s,[2 1])),'k^-');hold on;
        plot(te(:).',angle(permute(sHat(:),[2 1])),'x-');hold off;
        ylim([min(angle(s(:)))*0.9 max(angle(s(:)))*1.1]);
        title('Phase');
        
        subplot(2,3,3);
        plot(te(:).',real(permute(s,[2 1])),'k^-');hold on;
        plot(te(:).',imag(permute(s,[2 1])),'ks-');
        plot(te(:).',real(permute(sHat(:),[2 1])),'bx-.');
        plot(te(:).',imag(permute(sHat(:),[2 1])),'b*-.');hold off;
        ylim([min([real(s(:));imag(s(:))])*1.1,max([real(s(:));imag(s(:))])*1.1]);
        title('Real, Imaginary');
        
        subplot(2,3,4);
        plot(te(:).',real(permute(sHat(:),[2 1]))-real(permute(s(:),[2 1])),'rx-.');hold on;
        plot(te(:).',imag(permute(sHat(:),[2 1]))-imag(permute(s(:),[2 1])),'r*-.');hold off;
        title('Residual');
        
        ha = subplot(2,3,5); pos = get(ha,'Position'); un = get(ha,'Units'); delete(ha)
        uitable('Data',x(:),'Units',un,'Position',pos);
        
        subplot(2,3,6);
        plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm');
        text(0.5,0.5,sprintf('resnorm=%f',sum(err(:).^2)),'Units','normalized');
    end
    if length(DEBUG_resnormAll) <300
        xlim([0 300]);
    else
        xlim([length(DEBUG_resnormAll)-300 length(DEBUG_resnormAll)]);
    end
    drawnow;
end