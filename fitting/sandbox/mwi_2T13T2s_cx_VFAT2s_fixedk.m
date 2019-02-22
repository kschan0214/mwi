%% fitRes = mwi_2T13T2s_cx_VFAT2s(algoPara,imgPara)
%
% Input
% --------------
% algoPara.maxIter      : maximum iteration allows (default: 500)
% algoPara.fcnTol       : function tolerance (default: 1e-5)
% algoPara.isWeighted   : boolean cost function weighted by echo intensity (default: false)
% algoPara.weightMethod : if algoPara.isWeighted = true, then you may
%                         choose the weighting method (default: 'norm')(other: '1stEcho')
% algoPara.npulse       : no. of pulses to reach steady-state for EPG-X (default: 50)
% algoPara.numMagn      : no. of phase corrupted echoes (default: length(te))
%                         This will alter the type of fitting being used (magnitude/mixed/complex)
% algoPara.userDefine   : user defined x0,lb and ub (default: [])
% algoPara.isParallel   : parallel computing using parfor (default: false)
% algoPara.DEBUG        : debug mode (default: false)
% algoPara.verbose      : display feedback about the fitting (default: true)
%
% imgPara.img           : 5D image data, time in 4th dimension, T1w in 5th dim
% imgPara.mask          : signal mask
% imgPara.te            : echo times
% imgPara.fa            : flip angles
% imgPara.b1map         : B1 map
% imgPara.fieldmap      : field map
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
% Date last modified: 28 January 2019
%
%
function fitRes = mwi_2T13T2s_cx_VFAT2s_fixedk(algoPara,imgPara)
% make sure computeFiter in the same tool is used
% addpath('utils/');

disp('Myelin water imaing: VFA-T2* model');
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
npulse       = algoPara.npulse;
maxIter      = algoPara.maxIter;
fcnTol       = algoPara.fcnTol;
numMagn      = algoPara.numMagn;
isWeighted   = algoPara.isWeighted;
weightMethod = algoPara.weightMethod;
userDefine   = algoPara.userDefine;
isParallel   = algoPara.isParallel;
isInvivo     = algoPara.isInvivo;
model        = algoPara.model;

% check the user weighting method matches the available one or not
if strcmpi(weightMethod,'norm') && strcmpi(weightMethod,'1stEcho')
    disp('Do not support the input weighting method');
    weightMethod = 'norm';
end

% if DEBUG on then disables parallel computing
if DEBUG
    isParallel = false;
end

% capture all images related data
te    = imgPara.te;
tr    = imgPara.tr;
fa    = imgPara.fa;
data  = imgPara.img;
mask  = imgPara.mask;
b1map = imgPara.b1map;
fm    = imgPara.fieldmap;
pini  = imgPara.pini;

%
if isreal(data) && numMagn~=length(te)
    numMagn = length(te);
    disp('Input data is real, switch to magnitude data fitting');
end

%%%%%%%%%% display fitting message %%%%%%%%%%
if verbose
    disp('The following fitting parameters are used:');
    fprintf('Max. iterations = %i\n',maxIter);
    fprintf('Function tolerance = %e\n',fcnTol);
    
    switch model
        case 'epgx'
            fprintf('Model: EPG-X\n');
            fprintf('No. of pulses for EPG-X = %i\n',npulse);
            
        case 'epg'
            fprintf('Model: traditional EPG\n');
            fprintf('No. of pulses for EPG = %i\n',npulse);
            
        case 'standard'
            fprintf('Model: Standard\n');
            
    end
    
    if isWeighted
        disp('Weighted cost function: True');
        disp(['Weighting method: ' weightMethod]);
    else
        disp('Weighted cost function: False');
    end
    if isInvivo
        disp('Initial guesses for in vivo study');
    else
        disp('Initial guesses for ex vivo study');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny,nx,nz,~,nfa] = size(data);

% lsqnonlin setting
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*11,...
        'FunctionTolerance',fcnTol,'StepTolerance',1e-6);
  
% fdss = [1e-4,1e-4,1e-4,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8];  
% options.FiniteDifferenceStepSize = fdss;

% if DEBUG then display fitting message
if ~DEBUG
    options.Display = 'off';
end

if strcmpi(model,'epgx')
    % EPG-X has an extra exchange term
    if numMagn==numel(te)   % magnitude fitting has 11 estimates
        estimates = zeros(ny,nx,nz,10);
    else
        estimates = zeros(ny,nx,nz,10+nfa*2); % others have 12 estimates
    end
else
    if numMagn==numel(te)   % magnitude fitting has 10 estimates
        estimates = zeros(ny,nx,nz,10);
    else
        estimates = zeros(ny,nx,nz,11+nfa*2); % others have 12 estimates
    end
end

resnorm   = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        if verbose
            fprintf('Processing slice %i...\n',kz);
        end
%         numVoxelToBeFitted = length(find(mask(:,:,kz)==1));
%         fprintf('%i voxles need to be fitted...\n',numVoxelToBeFitted);
%         fprintf('Progress (%%): ');
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s       = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1      = b1map(ky,kx,kz);
                    db0     = squeeze(fm(ky,kx,kz,:));
                    pini0   = squeeze(pini(ky,kx,kz,:));
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,db0,pini0,npulse,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG);
                end
            end
%             if numVoxelToBeFitted > 0
%                 numVoxelFitted = length(find(mask(1:ky,:,kz)==1));
%                 percentFinish = floor(numVoxelFitted/numVoxelToBeFitted * 100);
%                 if percentFinish>1
%                   for j=0:log10(percentFinish-1)
%                       fprintf('\b'); % delete previous counter display
%                   end
%                   fprintf('%% %i ',percentFinish);
%                 end
%             end
            if mod(ky,5) == 0;
                fprintf('%i ', ky);
            end
        end
        fprintf('\n');
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
                    db0     = squeeze(fm(ky,kx,kz,:));
                    pini0   = squeeze(pini(ky,kx,kz,:));
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,db0,pini0,npulse,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG);
                end
            end
            if mod(ky,5) == 0;
                fprintf('%i ', ky);
            end
        end
        fprintf('\n');
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b1,db0,pini,npulse,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG)
% define initial guesses
% estimate rho0 of the first echo
[~,rho0] = DESPOT1(abs(s(:,1)),fa,tr,'b1',b1);
% estimate t1 from later echo
[t10,~] = DESPOT1(abs(s(:,end-3)),fa,tr,'b1',b1);

% in case rho0 and t10 go wrong then use defualt values
if rho0<0
    rho0=max(abs(s(:)));
end
if t10<0 || t10>2000e-3
    t10=1000e-3;
end

if isInvivo
    % common initial guesses for in vivo study
    Amw0   = 0.1*rho0; 	Amwlb   = 0;        Amwub   = 2*rho0;
    Aiw0   = 0.6*rho0; 	Aiwlb   = 0;        Aiwub   = 2*rho0;
    Aew0   = 0.3*rho0; 	Aewlb   = 0;        Aewub   = 2*rho0;
    t2smw0 = 10e-3;     t2smwlb = 1e-3;     t2smwub = 25e-3;
    t2siw0 = 64e-3; 	t2siwlb = 25e-3;    t2siwub = 500e-3;
    t2sew0 = 48e-3; 	t2sewlb = 25e-3;    t2sewub = 500e-3;
    t1s0   = 300e-3;  	t1slb   = 50e-3;  	t1sub   = 650e-3;
    t1l0   = t10;     	t1llb   = 500e-3; 	t1lub   = 5000e-3;
    kls0    = 2;       	klslb    = 0;      	klsub    = 10;       % exchange rate from long T1 to short T1

else
    % common initial guesses for ex vivo study
    Amw0   = 0.15*rho0; Amwlb   = 0;        Amwub   = 1*rho0;
    Aiw0   = 0.6*rho0;  Aiwlb   = 0;        Aiwub   = 1*rho0;
    Aew0   = 0.25*rho0; Aewlb   = 0;        Aewub   = 1*rho0;
    t2smw0 = 10e-3;     t2smwlb = 1e-3;     t2smwub = 20e-3;
    t2siw0 = 54e-3; 	t2siwlb = 20e-3;    t2siwub = 200e-3;
    t2sew0 = 38e-3; 	t2sewlb = 20e-3;    t2sewub = 200e-3;
    t1s0   = 100e-3;  	t1slb   = 50e-3;  	t1sub   = 400e-3;
    t1l0   = t10;     	t1llb   = 300e-3; 	t1lub   = 1500e-3;
    kls0   = 2;       	klslb   = 0;      	klsub   = 10;       % exchange rate from long T1 to short T1
    
end

    fmw0   = 5;   	fmwlb   = -75;  	fmwub   = 75;
    fiw0   = 0;    	fiwlb   = -25;      fiwub   = +25;
if numMagn==numel(te) % magnitude fitting
    
    % set initial guess and fitting boundaries
    x0 = double([Amw0 ,Aiw0 ,Aew0 ,t2smw0 ,t2siw0 ,t2sew0 ,t1s0 ,t1l0 ,fmw0 ,fiw0 ]);
    lb = double([Amwlb,Aiwlb,Aewlb,t2smwlb,t2siwlb,t2sewlb,t1slb,t1llb,fmwlb,fiwlb]);
    ub = double([Amwub,Aiwub,Aewub,t2smwub,t2siwub,t2sewub,t1sub,t1lub,fmwub,fiwub]);
    
    if strcmpi(model,'epgx')
        x0 = double([x0,kls0]);
        lb = double([lb,klslb]);
        ub = double([ub,klsub]);
    end
    
else    % other fittings
    totalField0 = db0(:).';  totalFieldlb   = db0(:).'-100;                 totalFieldub    = db0(:).'+100;
    pini0       = pini(:).'; pinilb         = ones(size(pini0))*(-2*pi);    piniub          = ones(size(pini0))*2*pi;
    
    if strcmpi(model,'epgx')
        % set initial guess and fitting boundaries
%         x0 = double([Amw0 ,Aiw0 ,Aew0 ,t2smw0 ,t2siw0 ,t2sew0 ,t1s0 ,t1l0 ,fmw0 ,fiw0  ,kls0,totalField0,pini0]);
%         lb = double([Amwlb,Aiwlb,Aewlb,t2smwlb,t2siwlb,t2sewlb,t1slb,t1llb,fmwlb,fiwlb,klslb,totalFieldlb,pinilb]);
%         ub = double([Amwub,Aiwub,Aewub,t2smwub,t2siwub,t2sewub,t1sub,t1lub,fmwub,fiwub,klsub,totalFieldub,piniub]);
        x0 = double([Amw0 ,Aiw0 ,Aew0 ,t2smw0 ,t2siw0 ,t2sew0 ,t1s0 ,t1l0 ,fmw0 ,fiw0  ,totalField0,pini0]);
        lb = double([Amwlb,Aiwlb,Aewlb,t2smwlb,t2siwlb,t2sewlb,t1slb,t1llb,fmwlb,fiwlb,totalFieldlb,pinilb]);
        ub = double([Amwub,Aiwub,Aewub,t2smwub,t2siwub,t2sewub,t1sub,t1lub,fmwub,fiwub,totalFieldub,piniub]);
    else
        % set initial guess and fitting boundaries
        x0 = double([Amw0 ,Aiw0 ,Aew0 ,t2smw0 ,t2siw0 ,t2sew0 ,t1s0 ,t1l0 ,fmw0 ,fiw0  ,totalField0,pini0]);
        lb = double([Amwlb,Aiwlb,Aewlb,t2smwlb,t2siwlb,t2sewlb,t1slb,t1llb,fmwlb,fiwlb,totalFieldlb,pinilb]);
        ub = double([Amwub,Aiwub,Aewub,t2smwub,t2siwub,t2sewub,t1sub,t1lub,fmwub,fiwub,totalFieldub,piniub]);

    end
    
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

% if DEBUG then create an array to store resnorm of all iterations
if DEBUG
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
    x0
end

% precompute EPG-X's transition matrix here for speed
RFphi0 = 50;      % initial RF phase, degrees
T3D_all = cell(length(fa),1);
if strcmpi(model,'epgx')
    phiCycle = RF_phase_cycle(npulse,RFphi0);
    for kfa=1:length(fa)
        T3D_all{kfa} = PrecomputeT(phiCycle,d2r(fa(kfa)*b1));
    end
elseif strcmpi(model,'epg')
    phiCycle = RF_phase_cycle(npulse,RFphi0);
    for kfa=1:length(fa)
        T3D_all{kfa} = precomputeT_epg(phiCycle,d2r(fa(kfa)*b1));
    end
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,npulse,RFphi0,T3D_all,numMagn,isWeighted,weightMethod,model,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,npulse,RFphi0,T3D_all,numMagn,isWeighted,weightMethod,model,DEBUG)
% capture all fitting parameters
Amw=x(1);   Aiw=x(2);   Aew=x(3);
t2smw=x(4); t2siw=x(5); t2sew=x(6);
t1s=x(7); t1l = x(8);
fmw=x(9);   fiw=x(10);  

few = 0; 
kls = 2;

if numMagn==numel(te) % magnitude fitting
           
%     if strcmpi(model,'epgx')
%         kls=x(10);
%     else
%         kls = 0;
%     end
    % no initial phase 
    pini=0;
    totalfield = 0;
else    % other fittings
    
    if strcmpi(model,'epgx')
%         kls=x(10);
        totalfield = x(11:11+length(fa)-1);
        pini       = x(11+length(fa):end);
    else
        kls=0;
        totalfield = x(10:10+length(fa)-1);
        pini       = x(10+length(fa):end);
    end
end

% simulate signal based on parameter input
% sHat = mwi_model_2T13T2scc_epgx(fa,te,tr,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1s,t1l,fmw,fiw,few,totalfield,pini,b1,kls,npulse,T3D_all);
sHat = mwi_model_2T13T2scc_epgx_moreoption(fa,te,tr,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1s,t1l,fmw,fiw,few,totalfield,pini,b1,kls,npulse,RFphi0,T3D_all,model);

% compute fitting residual
if isWeighted
    switch weightMethod
        case 'norm'
            % weights using echo intensity, as suggested in Nam's paper
%             w = abs(s(:))/norm(abs(s(:)));
            % in this way the sum of all weights is the same as no
            % weighting (which is sum(ones(size(s(:)))).)
%             w = numel(s) * abs(s(:))/sum(abs(s(:)));

            w = sqrt(abs(s)/sum(abs(s(:))));
%            w =  abs(s)/sum(abs(s(:)));
        case '1stEcho'
            % weights using the 1st echo intensity of each flip angle
            w = bsxfun(@rdivide,abs(s),abs(s(:,1)));
%             w = numel(s) * w(:)/sum(w(:));
            w = numel(s) * w/sum(w(:));
    end
    % compute the cost with weights
    err = computeFiter(s,sHat,numMagn,w);
else
    % compute the cost without weights
    err = computeFiter(s,sHat,numMagn);
end

% residual normalied by measured signal
% 20180316 TODO:maybe for weighted cost there should be another way to do this 
err = err./norm(abs(s));

% if DEBUG then plots current fitting result
if DEBUG
    global DEBUG_resnormAll
    figure(99);subplot(211);plot(te(:).',abs(permute(s,[2 1])),'^-');hold on;ylim([0-min(abs(s(:))),max(abs(s(:)))+10]);
    title('Magnitude');
    plot(te(:).',abs(permute(sHat,[2 1])),'x-.');plot(te(:).',(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'o-.');
    hold off;
    text(te(1)/3,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
    text(te(1)/3,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f,T1my=%f,T1l=%f,kmy=%f',...
        Amw,Aiw,Aew,t2smw,t2siw,t2sew,fmw,fiw,t1s,t1l,kls));
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
% copy input to output
imgPara2 = imgPara;
algoPara2 = algoPara;
isValid = true;

% check if the number of flip angles matches with the data's 5th dim
if length(imgPara.fa) ~= size(imgPara.img,5)
    isValid = false;
end
% check if the number of echo times matches with the data's 4th dim
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
    algoPara2.weightMethod = algoPara.weightMethod;
catch
    algoPara2.weightMethod = 'norm';
end
% check method for weighting
try
    algoPara2.isWeighted = algoPara.isWeighted;
catch
    algoPara2.isWeighted = false;
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
% check initial guesses
try
    imgPara2.isInvivo = imgPara.isInvivo;
catch
    imgPara2.isInvivo = true;
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
% check field map
try
    imgPara2.fieldmap = imgPara.fieldmap;
    disp('Field map input: True');
catch
    imgPara2.fieldmap = zeros(size(imgPara2.mask));
    disp('Field map input: False');
end
% check initial phase map
try
    imgPara2.pini = imgPara.pini;
    disp('Initial map input: True');
catch
    imgPara2.pini = zeros(size(imgPara2.mask));
    disp('Initial map input: False');
end

try
    algoPara2.model = algoPara.model;
    if ~strcmpi(algoPara2.model,'epgx') && ~strcmpi(algoPara2.model,'epg') && ~strcmpi(algoPara2.model,'standard')
        error('Your input model is not supported. Please choose either ''epgx'', ''epg'' or ''standard''. ');
    end
catch
    algoPara2.model = 'epgx';
end


end