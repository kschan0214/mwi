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
% Date created: 3 March 2019
% Date last modified: 
%
%
function fitRes = mwi_2T13T2s_cx_VFAT2s_HCFM(algoPara,imgPara)
disp('Myelin water imaing: VFA-T2* model with HCFM constraints');

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara]=CheckAndSetPara(algoPara,imgPara);

% get debug mode and verbose
DEBUG   = algoPara.DEBUG;

% capture all fitting settings
RFphase      = algoPara.RFphase;
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

% capture all images related data
te    = imgPara.te;
tr    = imgPara.tr;
fa    = imgPara.fa;
data  = imgPara.img;
mask  = imgPara.mask;
b1map = imgPara.b1map;
fm    = imgPara.fieldmap;
pini  = imgPara.pini;
intra = imgPara.intra;


[ny,nx,nz,~,nfa] = size(data);

% lsqnonlin setting
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*11,...
        'FunctionTolerance',fcnTol,'StepTolerance',1e-6);

if DEBUG
    % if DEBUG is on then disables parallel computing
    isParallel = false;
else
    % if not DEBUG then turn off fitting message
    options.Display = 'off';
end

if strcmpi(model,'epgx')
    % EPG-X has an extra exchange term
    if numMagn==numel(te)   % magnitude fitting has 11 estimates
        estimates = zeros(ny,nx,nz,7);
    else
        estimates = zeros(ny,nx,nz,8+nfa); % others have 12 estimates
    end
else
    if numMagn==numel(te)   % magnitude fitting has 10 estimates
        estimates = zeros(ny,nx,nz,7);
    else
        estimates = zeros(ny,nx,nz,7+nfa); % others have 12 estimates
    end
end

numMaskedVoxel = length(mask(mask==1));
numFittedVoxel = 0;
fprintf('%i voxles need to be fitted...\n',numMaskedVoxel);
fprintf('Progress (percent): ');
progress='';

resnorm   = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s       = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1      = b1map(ky,kx,kz);
                    db0     = squeeze(fm(ky,kx,kz,:));
                    pini0   = pini(ky,kx,kz);
                    intra0  = intra(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,intra0,db0,pini0,RFphase,npulse,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG);
                end
            end
            % display progress
            tmp = mask(ky,:,kz);
            numFittedVoxel = numFittedVoxel + length(tmp(tmp==1));
            for ii=1:length(progress); fprintf('\b'); end
            progress=sprintf('%0.1f', numFittedVoxel*100/numMaskedVoxel);
            fprintf(progress);
            
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
                    pini0   = pini(ky,kx,kz);
                    intra0  = intra(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,intra0,db0,pini0,RFphase,npulse,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG);
                end
            end
            % display progress
            tmp = mask(ky,:,kz);
            numFittedVoxel = numFittedVoxel + length(tmp(tmp==1));
            for ii=1:length(progress); fprintf('\b'); end
            progress=sprintf('%0.1f', numFittedVoxel*100/numMaskedVoxel);
            fprintf(progress);
        end
    end
end
fprintf('\n');

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b1,intra0,db0,pini,RFphase,npulse,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG)
% define initial guesses
% estimate rho0 of the first echo
[~,rho0] = DESPOT1(abs(s(:,1)),fa,tr,'b1',b1);

% in case rho0 and t10 go wrong then use defualt values
if rho0<0
    rho0=max(abs(s(:)));
end

if isInvivo
    % common initial guesses for in vivo study
    Amw0   = 0.1*rho0; 	Amwlb   = 0;        Amwub   = 2*rho0;
    Aiw0   = 0.6*rho0; 	Aiwlb   = 0;        Aiwub   = 2*rho0;
    Aew0   = 0.3*rho0; 	Aewlb   = 0;        Aewub   = 2*rho0;
    t2smw0 = 10;        t2smwlb = 1;        t2smwub = 25;
    t2fw0  = 64;        t2fwlb  = 25;       t2fwub  = 250;
    t1s0   = 118e-3;  	t1slb   = 50e-3;  	t1sub   = 650e-3;
    kls0    = 2;       	klslb    = 0;      	klsub    = 6;       % exchange rate from long T1 to short T1

else
    % common initial guesses for ex vivo study
    Amw0   = 0.15*rho0; Amwlb   = 0;        Amwub   = 1*rho0;
    Aiw0   = 0.6*rho0;  Aiwlb   = 0;        Aiwub   = 1*rho0;
    Aew0   = 0.25*rho0; Aewlb   = 0;        Aewub   = 1*rho0;
    t2smw0 = 10;        t2smwlb = 1;        t2smwub = 20;
    t2fw0  = 54;        t2fwlb  = 20;       t2fwub  = 200;
    t1s0   = 100e-3;  	t1slb   = 50e-3;  	t1sub   = 400e-3;
    kls0   = 2;       	klslb   = 0;      	klsub   = 6;       % exchange rate from long T1 to short T1
    
end
sin2theta0 = 0.5; sin2thetalb = 0; sin2thetaub = 1;

% set initial guess and fitting boundaries
x0 = double([Amw0 ,Aiw0 ,t2smw0 ,t2fw0 ,t1s0 ,sin2theta0 ]);
lb = double([Amwlb,Aiwlb,t2smwlb,t2fwlb,t1slb,sin2thetalb]);
ub = double([Amwub,Aiwub,t2smwub,t2fwub,t1sub,sin2thetaub]);

% optional exchange term
if strcmpi(model,'epgx')
    x0 = double([x0,kls0]);
    lb = double([lb,klslb]);
    ub = double([ub,klsub]);
end
% optional background fields and initial phase
if numMagn~=numel(te) % other fittings
      
    totalField0 = db0(:).';  totalFieldlb   = db0(:).'-100;	totalFieldub    = db0(:).'+100;
    pini0       = pini;      pinilb         = -2*pi;        piniub          = 2*pi;
    
    x0 = double([x0,totalField0 ,pini0]);
    lb = double([lb,totalFieldlb,pinilb]);
    ub = double([ub,totalFieldub,piniub]);
    
end
% optional extracellular water signal
if ~isnan(intra0)
    x0 = [x0,Aew0];
    lb = [lb,Aewlb];
    ub = [ub,Aewub];
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
phiCycle = RF_phase_cycle(npulse,RFphase);
if strcmpi(model,'epgx')
    for kfa=1:length(fa)
        T3D_all{kfa} = PrecomputeT(phiCycle,d2r(fa(kfa)*b1));
    end
    
elseif strcmpi(model,'epg')
    for kfa=1:length(fa)
        T3D_all{kfa} = precomputeT_epg(phiCycle,d2r(fa(kfa)*b1));
    end
    
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,intra0,npulse,RFphase,T3D_all,numMagn,isWeighted,weightMethod,model,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,intra,npulse,RFphi0,T3D_all,numMagn,isWeighted,weightMethod,model,DEBUG)
% capture all fitting parameters
Amw=x(1);   Aiw=x(2);   
t2smw=x(3); t2fw=x(4);
t1s=x(5);   
sin2theta = x(6);
if strcmpi(model,'epgx')
    kls=x(7);
else
    kls = 0;
end

if numMagn==numel(te) % magnitude fitting     
    % initial phase and background field = 0
    pini=0;
    freq_bkg = zeros(length(fa));
    
else    % other fittings      
    if strcmpi(model,'epgx')
        freq_bkg = x(8:8+length(fa)-1);
        pini       = x(8+length(fa));
    else
        freq_bkg = x(7:7+length(fa)-1);
        pini       = x(7+length(fa));
    end
end

if ~isnan(intra)
    Aew=Aiw*(1-intra)/intra;
else
    Aew=x(end);
end
%%%%%%%%%% 1. SPGR steady state %%%%%%%%%%
% simulate signal based on parameter input
%  [~,ssHat] = mwi_model_ssSPGR_2T1(fa*b1,tr,...
%                                     Amw,Aiw+Aew,...
%                                     t1s,t1l,...
%                                     [],kls);

rho_mw  = 0.43;
x_i     = -0.1;
x_a     = -0.1;
E       = 0;
b0      = 3;
gyro    = 42.57747892;
t1l     = 1.05;

g   = sqrt(Aiw/(Aiw+Amw/rho_mw));
if Aew == 0
    fvf = 1;
else
    fvf = (Aiw/Aew)/((Aiw/Aew)+g.^2);
end
% coefficients, Eq.[A15]
c1 = 1/4 - (3/2)*((g^2)/(1-g^2))*log(1/g);
%
fs = ((x_i/2)*(2/3-sin2theta) + (x_a/2)*(c1*sin2theta-1/3) + E) * gyro * b0;
% 
phiCycle = RF_phase_cycle(npulse,RFphi0);
t1x = [t1l, t1s]; 
t2x = [t2fw,t2smw]; % assuming T2* of iw has similar T2 of long T1 compartment
fx = Amw/(Aiw+Aew+Amw); % mwf
% fs = (freq_mw); % frequency difference between long and short T1 compartments
% fs = (fmw-(fiw*Aiw+few*Aew)/(Aiw+Aew)); % frequency difference between long and short T1 compartments
ssHat = zeros(2,length(fa));
for ii=1:length(fa) 
    % Compute RF spoling phase cycles
    % 2 pools, with exchange
    % start with steady-state signals, transverse magnitisation
    s1 = Signal_GRE_T1wMono((1-fx), fa(ii)*b1, t1x(1), tr);
    s2 = Signal_GRE_T1wMono(fx, fa(ii)*b1, t1x(2), tr);
    
    % EPG-X core
    z1 = s1/sind(fa(ii)*b1);    % long T1 compartment, longitudinal magnitisation (z-state)
    z2 = s2/sind(fa(ii)*b1);    % short T1 compartment, longitudinal magnitisation (z-state)
    tmp = EPGX_GRE_BMsplit_PrecomputedT(T3D_all{ii},phiCycle,tr,t1x,t2x,fx,kls,'delta',fs,'kmax',10,'ss',[z1,z2]);
    ssHat(1,ii) = (tmp{2}(end))*(Aiw+Aew+Amw);
    ssHat(2,ii) = (tmp{1}(end))*(Aiw+Aew+Amw);
     
end                                           

%%%%%%%%%% 2. T2* decay %%%%%%%%%%
param.b0        = 3;
param.pini      = pini;
param.g         = g;
param.fvf       = fvf;

sHat = zeros(size(s));
for kfa = 1:length(fa)
    param.freq_bkg  = freq_bkg(kfa);
    
    sHat(kfa,:) = mwi_model_ssSPGR_3T2scc_HCFM(te,ssHat(1,kfa),ssHat(2,kfa)*(Aiw/(Aiw+Aew)),ssHat(2,kfa)*(Aew/(Aiw+Aew)),...
                                     t2smw,t2fw,...
                                     sin2theta,param);
end

% compute fitting residual
if isWeighted
    switch weightMethod
        case 'norm'
            w = sqrt(abs(s)/sum(abs(s(:))));
            
        case '1stEcho'
            % weights using the 1st echo intensity of each flip angle
            w = bsxfun(@rdivide,abs(s),abs(s(:,1)));
%             w = numel(s) * w(:)/sum(w(:));
            w = numel(s) * w/sum(w(:));
    end
    
else
    % compute the cost without weights
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

%% check and set default
function [algoPara2,imgPara2]=CheckAndSetPara(algoPara,imgPara)
% copy input to output
imgPara2 = imgPara;
algoPara2 = algoPara;

%%%%%%%%%% 1. check algorithm parameters %%%%%%%%%%
% check debug
try algoPara2.DEBUG = algoPara.DEBUG;                   catch; algoPara2.DEBUG = false; end
% check parallel computing 
try algoPara2.isParallel = algoPara.isParallel;         catch; algoPara2.isParallel = false; end
% check RF phase for EPG/EPGX
try algoPara2.RFphase = algoPara.RFphae;                catch; algoPara2.RFphae = 50; end
% check number of pulse for EPG/EPGX
try algoPara2.npulse = algoPara.npulse;                 catch; algoPara2.npulse = 200; end
% check maximum iterations allowed
try algoPara2.maxIter = algoPara.maxIter;               catch; algoPara2.maxIter = 500; end
% check function tolerance
try algoPara2.fcnTol = algoPara.fcnTol;                 catch; algoPara2.fcnTol = 1e-5; end
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
    imgPara2.pini = zeros([size(imgPara2.mask) size(imgPara2.fa)]);
    disp('Initial map input: False');
end
% check intra-neurite volume fraction map
try
    imgPara2.intra = imgPara.intra;
    disp('Intra-neurite volume fraction map input: True');
catch
    imgPara2.intra = ones(size(imgPara2.mask))*nan;
    disp('Intra-neurite volume fraction map input: False');
end

% check input image is real or not
if isreal(imgPara2.data) && algoPara2.numMagn~=length(imgPara2.te)
    algoPara2.numMagn = length(imgPara2.tete);
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
fprintf('Function tolerance = %e\n',algoPara2.fcnTol);
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
    disp('Only ''norm'' or ''1stEcho'' are supported. Using ''norm '' instead');
    algoPara2.weightMethod = 'norm';
end

end

function Debug_display(s,sHat,err,te,fa,x,numMagn)
global DEBUG_resnormAll
DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    figure(99);

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
        ylim([min(angle(s(:)))*0.9 max(angle(s(:)))*1.1]);
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
end