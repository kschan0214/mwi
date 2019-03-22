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
function fitRes = mwi_2T13T2s_cx_VFAT2s_2steps_varpro(algoPara,imgPara)

disp('Myelin water imaing: VFA-T2* model');
% check validity of the algorithm parameters and image parameters
[algoPara,imgPara]=CheckAndSetDefault(algoPara,imgPara);

% get debug mode and verbose
DEBUG   = algoPara.DEBUG;

% capture all fitting settings
RFphase         = algoPara.RFphase;
npulse          = algoPara.npulse;
maxIter         = algoPara.maxIter;
fcnTol          = algoPara.fcnTol;
stepTol         = algoPara.stepTol;
numMagn         = algoPara.numMagn;
isWeighted      = algoPara.isWeighted;
weightMethod    = algoPara.weightMethod;
userDefine      = algoPara.userDefine;
isParallel      = algoPara.isParallel;
isInvivo        = algoPara.isInvivo;
model           = algoPara.model;
numEst_step1	= algoPara.numEst_step1;
numEst_step2  	= algoPara.numEst_step2;


% capture all images related data
te    = imgPara.te;
tr    = imgPara.tr;
fa    = imgPara.fa;
data  = imgPara.img;
mask  = imgPara.mask;
b1map = imgPara.b1map;
fm    = imgPara.fieldmap;
pini  = imgPara.pini;

[ny,nx,nz,~,~] = size(data);

% set fitting options
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*(numEst_step1+numEst_step2),...
    'StepTolerance',stepTol,'FunctionTolerance',fcnTol);

if DEBUG
    % if DEBUG is on then disables parallel computing
    isParallel = false;
else
    options.Display = 'off';
end

% create empty array for fitting results
estimates_step1 = zeros(ny,nx,nz,numEst_step1);
estimates_step2 = zeros(ny,nx,nz,numEst_step2);

numMaskedVoxel = length(mask(mask==1));
numFittedVoxel = 0;
fprintf('%i voxel(s) to be fitted...\n',numMaskedVoxel);
progress='0 %%';
fprintf('Progress: ');
fprintf(progress);

resnorm_step1 = zeros(ny,nx,nz);
resnorm_step2 = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s       = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1      = b1map(ky,kx,kz);
                    db0     = squeeze(fm(ky,kx,kz,:));
                    pini0   = squeeze(pini(ky,kx,kz,:));
                    [estimates_step1(ky,kx,kz,:),resnorm_step1(ky,kx,kz)] = FitModel_t2s(s,te,db0,numMagn,isWeighted,weightMethod,isInvivo,userDefine,options,DEBUG);
%                     [estimates_step2(ky,kx,kz,:),resnorm_step2(ky,kx,kz)] = FitModel_vfa(s,fa,te,tr,b1,pini0,squeeze(estimates_step1(ky,kx,kz,:)),npulse,RFphase,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG);
                end
            end
            % display progress
            [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,mask(ky,:,kz));
            
        end
    end
else
    for kz=1:nz
        for ky=1:ny
            for kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1 = b1map(ky,kx,kz);
                    db0     = squeeze(fm(ky,kx,kz,:));
                    pini0   = squeeze(pini(ky,kx,kz,:));
                    [estimates_step1(ky,kx,kz,:),resnorm_step1(ky,kx,kz)] = FitModel_t2s(s,te,db0,numMagn,isWeighted,weightMethod,isInvivo,userDefine,options,DEBUG);
%                     [estimates_step2(ky,kx,kz,:),resnorm_step2(ky,kx,kz)] = FitModel_vfa(s,fa,te,tr,b1,pini0,squeeze(estimates_step1(ky,kx,kz,:)),npulse,RFphase,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG);
                end
            end
            % display progress
            [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,mask(ky,:,kz));
            
        end
    end
end

fitRes.estimates_step1 = estimates_step1;
fitRes.resnorm_step1   = resnorm_step1;
% fitRes.estimates_step2 = estimates_step2;
% fitRes.resnorm_step2   = resnorm_step2;
% fitRes.estimates       = cat(4,estimates_step2(:,:,:,1:3),...
%                                estimates_step1(:,:,:,1:3),...
%                                estimates_step2(:,:,:,4:5),...
%                                estimates_step1(:,:,:,4:5),...
%                                estimates_step2(:,:,:,6:end));

end

%% Setup lsqnonlin and fit with the MC-T2s model
function [x,res] = FitModel_t2s(s,te,db0,numMagn,isWeighted,weightMethod,isInvivo,userDefine,options,DEBUG)

if isInvivo
    % common initial guesses for in vivo study
    t2smw0 = 10e-3;            	t2smwlb = 1e-3;               	t2smwub = 25e-3;
    t2siw0 = 64e-3;            	t2siwlb = 25e-3;                t2siwub = 200e-3;
    t2sew0 = 48e-3;           	t2sewlb = 25e-3;                t2sewub = 200e-3;
    
else
    % common initial guesses for ex vivo study
    t2smw0 = 10e-3;                 t2smwlb = 1e-3;                	t2smwub = 20e-3;
    t2siw0 = 54e-3;                 t2siwlb = 20e-3;                t2siwub = 200e-3;
    t2sew0 = 38e-3;                 t2sewlb = 20e-3;                t2sewub = 200e-3;
    
end

% set initial guess and fitting boundaries
fmw0   = 5;   	fmwlb   = -75;  	fmwub   = +75;
fiw0   = -2;    	fiwlb   = -25;      fiwub   = +25;

x0 = double([t2smw0 ,t2siw0 ,t2sew0 ,fmw0 ,fiw0 ]);
lb = double([t2smwlb,t2siwlb,t2sewlb,fmwlb,fiwlb]);
ub = double([t2smwub,t2siwub,t2sewub,fmwub,fiwub]);

if numMagn~=numel(te) % other fittings
       
    % total field in Hz, initial phase in radian
    totalField0 = db0(:).';  totalFieldlb   = db0(:).'-25;                 totalFieldub    = db0(:).'+25;
    
    % set initial guess and fitting boundaries
    x0 = double([x0 ,totalField0]);
    lb = double([lb,totalFieldlb]);
    ub = double([ub,totalFieldub]);
    
end

% set user defined initial guess and fitting bounds here
if ~isempty(userDefine.x0)
    x0(~isnan(userDefine.x0)) = userDefine.x0(~isnan(userDefine.x0));
end
if ~isempty(userDefine.lb)
    lb(~isnan(userDefine.lb)) = userDefine.lb(~isnan(userDefine.lb));
end
if ~isempty(userDefine.ub)
    ub(~isnan(userDefine.ub)) = userDefine.ub(~isnan(userDefine.ub));
end

% if DEBUG then create an array to store resnorm of all iterations
if DEBUG
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
    x0
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc_t2s(y,s,te,numMagn,isWeighted,weightMethod,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc_t2s(x,s,te,numMagn,isWeighted,weightMethod,DEBUG)
nfa = size(s,1);
% capture all fitting parameters
t2smw=x(1);   t2siw=x(2);   t2sew=x(3);
fmw=x(4);     fiw=x(5); 
few = 0;  

if numMagn==numel(te) % magnitude fitting
    % no initial phase 
    totalfield  = zeros(1,nfa);
    
else    % other fittings
    totalfield = x(6:end);
end

% phi matrix is time-dependent
sHat = zeros(size(s));
phi_mat = [-1/t2smw+1i*2*pi*fmw, -1/t2siw+1i*2*pi*fiw, -1/t2sew+1i*2*pi*few];
phi_mat = exp(te(:) * phi_mat);
for kfa = 1:nfa
    lambda_mat = diag(exp(1i*2*pi*totalfield(kfa) * te(:)));
    psi_mat = lambda_mat * phi_mat;
    sHat(kfa,:) = psi_mat*pinv(psi_mat)*s(kfa,:).';
%     lambda_mat_inv = diag(exp(-1i*2*pi*totalfield(kfa) * te(:)));
%     sHat(kfa,:) = psi_mat*((phi_mat'*phi_mat)\phi_mat'*lambda_mat_inv)*s(kfa,:).';
end

% compute fitting residual
if isWeighted
    switch weightMethod
        case 'norm'
            % weights using echo intensity, as suggested in Nam's paper
            w = sqrt(abs(s)/sum(abs(s(:))));

        case '1stEcho'
            % weights using the 1st echo intensity of each flip angle
            w = bsxfun(@rdivide,abs(s),abs(s(:,1)));
%             w = numel(s) * w(:)/sum(w(:));
            w = numel(s) * w/sum(w(:));
    end
else
    % no weighting = same weighting
    w = sqrt(ones(size(s))/numel(s));
end
% compute the cost with weights
err = computeFiter(s,sHat,numMagn,w);
    
% residual normalied by measured signal
% 20180316 TODO:maybe for weighted cost there should be another way to do this 
% err = err./norm(abs(s));

% if DEBUG then plots current fitting result
if DEBUG
    Debug_display(s,sHat,err,te,fa,x,numMagn)
end

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel_vfa(s,fa,te,tr,b1,pini,estimates,npulse,RFphi0,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG)

% estimate proton density of total water using 1st echo
[~,rho0] = DESPOT1(abs(s(:,1)),fa,tr,'b1',b1);
% estimate long t1 from later echoes
[t10,~] = DESPOT1(abs(s(:,end-3)),fa,tr,'b1',b1);

% in case rho0 and t10 go wrong then use defualt values
if rho0<0
    rho0=max(abs(s(:)));
end
if t10<0 || t10>3000e-3
    t10=1000e-3;
end

Amw0   = 0.1*rho0; 	Amwlb   = 0;        Amwub   = 2*rho0;
Aiw0   = 0.6*rho0; 	Aiwlb   = 0;        Aiwub   = 2*rho0;
Aew0   = 0.4*rho0; 	Aewlb   = 0;        Aewub   = 2*rho0;
if isInvivo
    % common initial guesses for in vivo study
    t1s0   = 118e-3;  	t1slb   = 50e-3;  	t1sub   = 650e-3;
    t1l0   = t10;     	t1llb   = 500e-3; 	t1lub   = 3000e-3;
%     t1s0   = t1(1);         t1slb   = 50e-3;  	t1sub   = 650e-3;
%     t1l0   = mean(t1(2:3)); t1llb   = 500e-3; 	t1lub   = 3000e-3;
%     t1s0   = 500e-3;         t1slb   = 50e-3;  	t1sub   = 650e-3;
%     t1l0   = 3; t1llb   = 500e-3; 	t1lub   = 3000e-3;
    kls0    = 2;       	klslb    = 0;      	klsub    = 6;       % exchange rate from long T1 to short T1

else
    % common initial guesses for ex vivo study
    t1s0   = 100e-3;  	t1slb   = 50e-3;  	t1sub   = 400e-3;
    t1l0   = t10;     	t1llb   = 300e-3; 	t1lub   = 2000e-3;
    kls0   = 2;       	klslb   = 0;      	klsub   = 6;       % exchange rate from long T1 to short T1
    
end

% set initial guess and fitting boundaries
x0 = double([Amw0 ,Aiw0 ,Aew0 ,t1s0 ,t1l0]);
lb = double([Amwlb,Aiwlb,Aewlb,t1slb,t1llb]);
ub = double([Amwub,Aiwub,Aewub,t1sub,t1lub]);

if strcmpi(model,'epgx')
    x0 = double([x0,kls0]);
    lb = double([lb,klslb]);
    ub = double([ub,klsub]);
end

if numMagn~=numel(te) % other fittings
   
    pini0       = pini;   	pinilb         = -2*pi;         piniub          = 2*pi;

    % set initial guess and fitting boundaries
    x0 = double([x0,pini0]);
    lb = double([lb,pinilb]);
    ub = double([ub,piniub]);
    
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

% if DEBUG then create an array to store resnorm of all iterations
if DEBUG
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
    x0
end

% precompute EPG-X's transition matrix here for speed
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
[x,res] = lsqnonlin(@(y)CostFunc_vfa(y,s,fa,te,tr,b1,estimates,npulse,RFphi0,T3D_all,numMagn,isWeighted,weightMethod,model,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc_vfa(x,s,fa,te,tr,b1,estimates,npulse,RFphi0,T3D_all,numMagn,isWeighted,weightMethod,model,DEBUG)
% fixed parameters
t2smw=estimates(1); t2siw=estimates(2); t2sew=estimates(3);
fmw=estimates(4);   fiw=estimates(5); 
if length(estimates) > 5
    totalfield = estimates(6:end);
else
    totalfield = zeros(1,length(fa));
end
few = 0;

% fitting parameters
Amw=x(1);   Aiw=x(2);   Aew=x(3);
t1s=x(4);   t1l=x(5);

if strcmpi(model,'epgx')
    kls=x(6);
else
    kls = 0;
end
 
if numMagn==numel(te) % magnitude fitting
    % no initial phase 
    pini=0;
    
else    % other fittings  
    if strcmpi(model,'epgx')
        pini       = x(7);
    else
        pini       = x(6);
    end
end

% simulate signal based on parameter input
sHat = mwi_model_2T13T2scc_epgx_moreoption(fa,te,tr,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1s,t1l,fmw,fiw,few,totalfield,pini,b1,kls,npulse,RFphi0,T3D_all,model);

% compute fitting residual
if isWeighted
    switch weightMethod
        case 'norm'
            % weights using echo intensity, as suggested in Nam's paper
            w = sqrt(abs(s)/sum(abs(s(:))));

        case '1stEcho'
            % weights using the 1st echo intensity of each flip angle
            w = bsxfun(@rdivide,abs(s),abs(s(:,1)));
%             w = numel(s) * w(:)/sum(w(:));
            w = numel(s) * w/sum(w(:));
    end
else
    % no weighting = same weighting
    w = sqrt(ones(size(s))/numel(s));
end
% compute the cost with weights
err = computeFiter(s,sHat,numMagn,w);

% residual normalied by measured signal
% 20180316 TODO:maybe for weighted cost there should be another way to do this 
err = err./norm(abs(s));

% if DEBUG then plots current fitting result
if DEBUG
    Debug_display(s,sHat,err,te,fa,x,numMagn)
end

end

% check and set default
function [algoPara2,imgPara2]=CheckAndSetDefault(algoPara,imgPara)
% copy input to output
imgPara2 = imgPara;
algoPara2 = algoPara;

%%%%%%%%%% 1. check algorithm parameters %%%%%%%%%%
% check debug
try algoPara2.DEBUG = algoPara.DEBUG;                   catch; algoPara2.DEBUG = false; end
% check parallel computing 
try algoPara2.isParallel = algoPara.isParallel;         catch; algoPara2.isParallel = false; end
% check RF phase for EPG/EPGX
try algoPara2.RFphase = algoPara.RFphae;                catch; algoPara2.RFphase = 50; end
% check number of pulse for EPG/EPGX
try algoPara2.npulse = algoPara.npulse;                 catch; algoPara2.npulse = 200; end
% check maximum iterations allowed
try algoPara2.maxIter = algoPara.maxIter;               catch; algoPara2.maxIter = 500; end
% check function tolerance
try algoPara2.fcnTol = algoPara.fcnTol;                 catch; algoPara2.fcnTol = 1e-5; end
% check step tolerance
try algoPara2.stepTol = algoPara.stepTol;               catch; algoPara2.stepTol = 1e-6; end
% check weighted sum of cost function
try algoPara2.weightMethod = algoPara.weightMethod;     catch; algoPara2.weightMethod = 'norm'; end
% check method for weighting
try algoPara2.isWeighted = algoPara.isWeighted;         catch; algoPara2.isWeighted = false; end
% check # of phase-corrupted echoes
try algoPara2.numMagn = algoPara.numMagn;               catch; algoPara2.numMagn = numel(imgPara.te); end
% check initial guesses
try algoPara2.isInvivo = algoPara2.isInvivo;            catch; algoPara2.isInvivo = true; end
% check model
try algoPara2.model = algoPara.model;                   catch; algoPara2.model = 'epgx'; end
% check user bounds and initial guesses
try algoPara2.userDefine.x0 = algoPara.userDefine.x0;   catch; algoPara2.userDefine.x0 = [];end
try algoPara2.userDefine.lb = algoPara.userDefine.lb;   catch; algoPara2.userDefine.lb = [];end
try algoPara2.userDefine.ub = algoPara.userDefine.ub;   catch; algoPara2.userDefine.ub = [];end

%%%%%%%%%% 2. check data integrity %%%%%%%%%%
% check if the number of flip angles matches with the data's 5th dim
if length(imgPara.fa) ~= size(imgPara.img,5)
    error('The length of flip angle does not match with the last dimension of the image.');
end
% check if the number of echo times matches with the data's 4th dim
if length(imgPara.te) ~= size(imgPara.img,4)
    error('The length of TE does not match with the 4th dimension of the image.');
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
    imgPara2.b1map = ones(size(imgPara2.mask));
    disp('B1 input: False');
end
% check field map
try
    imgPara2.fieldmap = imgPara.fieldmap;
    disp('Field map input: True');
catch
    imgPara2.fieldmap = zeros([size(imgPara2.mask) size(imgPara2.fa)]);
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

% check input image is real or not
if isreal(imgPara2.img) && algoPara2.numMagn~=length(imgPara2.te)
    algoPara2.numMagn = length(imgPara2.te);
    warning('Input data is real, switched to magnitude data fitting');
end

%%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
disp('The following fitting parameters are used:');
switch algoPara2.model
    case 'epgx'
        fprintf('Model: EPG-X\n');
        fprintf('No. of pulses for EPG-X = %i\n',algoPara2.npulse);
        fprintf('RF phase = %i\n',algoPara2.RFphase);

    case 'epg'
        fprintf('Model: EPG with no exchange\n');
        fprintf('No. of pulses for EPG = %i\n',algoPara2.npulse);
        fprintf('RF phase = %i\n',algoPara2.RFphase);

    case 'standard'
        fprintf('Model: Standard\n');
        
    otherwise
        error('Your input model is not supported. Please choose either ''epgx'', ''epg'' or ''standard''. ');
        
end

disp('Fitting options:');
fprintf('Max. iterations = %i\n',algoPara2.maxIter);
fprintf('Function tolerance = %.2e\n',algoPara2.fcnTol);
fprintf('Step tolerance = %.2e\n',algoPara2.stepTol);
if algoPara2.isInvivo
    disp('Initial guesses for in vivo study');
else
    disp('Initial guesses for ex vivo study');
end
% type of fitting
if algoPara2.numMagn==0
    disp('Fitting complex model with complex data');
elseif algoPara2.numMagn==numel(imgPara2.te)
    disp('Fitting complex model with magnitude data');
else
    fprintf('Fitting complex model with %i magnitude data and %i complex data\n',algoPara2.numMagn,numel(imgPara2.te)-algoPara2.numMagn);
end
% initial guesses and fitting bounds
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

disp('Cost function options:');
if algoPara2.isWeighted
    disp('Weighted cost function: True');
    disp(['Weighting method: ' algoPara2.weightMethod]);
else
    disp('Weighted cost function: False');
end
% check the user weighting method matches the available one or not
if strcmpi(algoPara2.weightMethod,'norm') && strcmpi(algoPara2.weightMethod,'1stEcho')
    disp(['Do not support the input weighting method: ' algoPara2.weightMethod]);
    disp('Only ''norm'' and ''1stEcho'' are supported. Using ''norm '' instead');
    algoPara2.weightMethod = 'norm';
end

% determine the number of estimates
% determine the number of estimates
numEst_step1 = 5; % basic setting has 5 estimates
if algoPara2.numMagn~=numel(imgPara2.te)
    numEst_step1 = numEst_step1 + length(imgPara2.fa); % total field and inital phase
end
numEst_step2 = 5; % basic setting has 3 estimates
if algoPara2.numMagn~=numel(imgPara2.te)
    numEst_step2 = numEst_step2 + 1; % total field and inital phase
end
if strcmpi(algoPara2.model,'epgx')
    numEst_step2 = numEst_step2 + 1; % total field and inital phase
end
algoPara2.numEst_step1 = numEst_step1;
algoPara2.numEst_step2 = numEst_step2;

end

%% debug display function
function Debug_display(s,sHat,err,te,fa,x,numMagn)
global DEBUG_resnormAll
DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    figure(97);

    if numMagn==numel(te)
        subplot(2,2,1);
        plot(te(:).',abs(permute(s,[2 1])),'^-');hold on;
        plot(te(:).',abs(permute(sHat,[2 1])),'x-.');hold off;
        ylim([min(abs(s(:)))*0.9,max(abs(s(:)))*1.1]); 
        for kfa = 1:length(fa); text(te(1)/3,abs(s(kfa,1)),['FA ' num2str(fa(kfa))]); end
        title('Magn.');
        
        subplot(2,2,2); 
        plot(te(:).',(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'o-.'); 
        title('residual');
        
        ha = subplot(2,2,3); pos = get(ha,'Position'); un = get(ha,'Units'); delete(ha)
        uitable('Data',x(:),'Units',un,'Position',pos);
        
        subplot(2,2,4);
        plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm')
        text(0.5,0.5,sprintf('resnorm=%e',sum(err(:).^2)),'Units','normalized');
    else
        subplot(2,3,1);
        plot(te(:).',abs(permute(s,[2 1])),'^-');hold on;
        plot(te(:).',abs(permute(sHat,[2 1])),'x-.');hold off;
        ylim([min(abs(s(:)))*0.9,max(abs(s(:)))*1.1]);
        for kfa = 1:length(fa); text(te(1)/3,abs(s(kfa,1)),['FA ' num2str(fa(kfa))]); end
        title('Magn.');
        
        subplot(2,3,2);
        plot(te(:).',angle(permute(s,[2 1])),'^-');hold on;
        plot(te(:).',angle(permute(sHat,[2 1])),'x-');hold off;
        ylim([min(angle(s(:))) max(angle(s(:)))*1.1]);
        title('Phase');
        
        subplot(2,3,3);
        plot(te(:).',real(permute(s,[2 1])),'^-');hold on;
        plot(te(:).',imag(permute(s,[2 1])),'s-');
        plot(te(:).',real(permute(sHat,[2 1])),'x-.');
        plot(te(:).',imag(permute(sHat,[2 1])),'*-.');hold off;
        ylim([min([real(s(:));imag(s(:))]),max([real(s(:));imag(s(:))])*1.1]);
        title('Real, Imaginary');
        
        subplot(2,3,4);
        plot(te(:).',real(permute(sHat-s,[2 1])),'x-.');hold on;
        plot(te(:).',imag(permute(sHat-s,[2 1])),'*-.');hold off;
        title('Residual');
        
        ha = subplot(2,3,5); pos = get(ha,'Position'); un = get(ha,'Units'); delete(ha)
        uitable('Data',x(:),'Units',un,'Position',pos);
        
        subplot(2,3,6);
        plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm');
        text(0.5,0.5,sprintf('resnorm=%e',sum(err(:).^2)),'Units','normalized');
    end
    if length(DEBUG_resnormAll) <100
        xlim([0 100]);
    else
        xlim([length(DEBUG_resnormAll)-100 length(DEBUG_resnormAll)]);
    end
    drawnow;
end

%% progress display
function [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,mask)
% number of non zeros element in the current mask
numFittedVoxel = numFittedVoxel + nnz(mask);
% delete previous progress
for ii=1:length(progress)-1; fprintf('\b'); end
% display current progress
progress=sprintf('%d %%%%', floor(numFittedVoxel*100/numMaskedVoxel));
fprintf(progress);

end