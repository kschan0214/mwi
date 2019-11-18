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
% Date last modified: 08 march 2018 
%
%
function fitRes = mwi_2T13T2s_cx_testing(algoPara,imgPara)
disp('Myelin water imaing: VFA-T2* model');
% check validity of the algorithm parameters and image parameters
[algoPara,imgPara]=CheckAndSetPara(algoPara,imgPara);

% get debug mode and verbose
DEBUG   = algoPara.DEBUG;

% capture all fitting settings
npulse       = algoPara.npulse;
maxIter      = algoPara.maxIter;
fcnTol       = algoPara.fcnTol;
stepTol      = algoPara.stepTol;
numMagn      = algoPara.numMagn;
isWeighted   = algoPara.isWeighted;
weightMethod = algoPara.weightMethod;
userDefine   = algoPara.userDefine;
isParallel   = algoPara.isParallel;
isInvivo     = algoPara.isInvivo;
numEst       = algoPara.numEst;
isNormCost   = algoPara.isNormCost;

% check the user weighting method matches the available one or not
if strcmpi(weightMethod,'norm') && strcmpi(weightMethod,'1stEcho')
    disp('Do not support the input weighting method');
    weightMethod = 'norm';
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

[ny,nx,nz,~,nfa] = size(data);

% lsqnonlin setting
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*numEst,...
        'FunctionTolerance',fcnTol,'StepTolerance',stepTol);

if DEBUG
    % if DEBUG is on then disables parallel computing
    isParallel = false;
else
    options.Display = 'off';
end

estimates = zeros(ny,nx,nz,numEst);

numMaskedVoxel = length(mask(mask==1));
numFittedVoxel = 0;
fprintf('%i voxel(s) to be fitted...\n',numMaskedVoxel);
progress='0 %%';
fprintf('Progress: ');
fprintf(progress);

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
                    pini0   = squeeze(pini(ky,kx,kz,:));
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,db0,pini0,npulse,numMagn,isWeighted,weightMethod,isNormCost,isInvivo,userDefine,options,DEBUG);
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
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,db0,pini0,npulse,numMagn,isWeighted,weightMethod,isNormCost,isInvivo,userDefine,options,DEBUG);
                end
            end
            % display progress
            [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,mask(ky,:,kz));
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b1,db0,pini,npulse,numMagn,isWeighted,weightMethod,isNormCost,isInvivo,userDefine,options,DEBUG)
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

%%%%%%%%%% Step 1: determine and set initial guesses  %%%%%%%%%%
if isInvivo
    % range for in vivo
    r2smw0 = 100;	r2smwlb = 40;	r2smwub = 300;
    r2siw0 = 16;	r2siwlb = 6;	r2siwub = 40;
    r2sew0 = 21;	r2sewlb = 6;	r2sewub = 40;
    mwf = 0.1;
    iwf = 0.6;
    t1s0   = 225e-3;  	t1slb   = 50e-3;  	t1sub   = 650e-3;
    t1l0   = 1.07;     	t1llb   = 500e-3; 	t1lub   = 3000e-3;
    kls0    = 2;       	klslb    = 0;      	klsub    = 6;       % exchange rate from long T1 to short T1
else
    % range for ex vivo
    r2smw0 = 150;	r2smwlb = 40;	r2smwub = 300;
    r2siw0 = 20;	r2siwlb = 6;	r2siwub = 40;
    r2sew0 = 30;	r2sewlb = 6;	r2sewub = 40;
    mwf = 0.15;
    iwf = 0.6;
    t1s0   = 225e-3;  	t1slb   = 50e-3;  	t1sub   = 400e-3;
    t1l0   = 1.07;     	t1llb   = 300e-3; 	t1lub   = 1500e-3;
    kls0   = 2;       	klslb   = 0;      	klsub   = 6;       % exchange rate from long T1 to short T1
end



Amw0   = mwf*rho0; 	Amwlb   = 0;            Amwub   = 1.5*rho0;
Aiw0   = iwf*rho0; 	Aiwlb   = 0;            Aiwub   = 1.5*rho0;
Aew0   = (1-mwf-iwf)*rho0; 	Aewlb   = 0;  	Aewub   = 1.5*rho0;

fmw0   = 5;   	fmwlb   = 5-75;  	fmwub   = 5+75;
fiw0   = 0;    	fiwlb   = -25;      fiwub   = +25;

% x0 = double([Amw0 ,Aiw0 ,Aew0 ,r2smw0 ,r2siw0 ,r2sew0 ,t1s0 ,t1l0 ,fmw0 ,fiw0,kls0]);
% lb = double([Amwlb,Aiwlb,Aewlb,r2smwlb,r2siwlb,r2sewlb,t1slb,t1llb,fmwlb,fiwlb,klslb]);
% ub = double([Amwub,Aiwub,Aewub,r2smwub,r2siwub,r2sewub,t1sub,t1lub,fmwub,fiwub,klsub]);
% x0 = double([Amw0 ,Aiw0 ,Aew0 ,r2smw0 ,r2siw0 ,r2sew0 ,t1s0 ,t1l0 ,fmw0 ,fiw0]);
% lb = double([Amwlb,Aiwlb,Aewlb,r2smwlb,r2siwlb,r2sewlb,t1slb,t1llb,fmwlb,fiwlb]);
% ub = double([Amwub,Aiwub,Aewub,r2smwub,r2siwub,r2sewub,t1sub,t1lub,fmwub,fiwub]);
x0 = double([Amw0 ,Aiw0 ,Aew0 ,r2smw0 ,r2siw0 ,r2sew0 ,fmw0 ,fiw0,0,0]);
lb = double([Amwlb,Aiwlb,Aewlb,r2smwlb,r2siwlb,r2sewlb,fmwlb,fiwlb,0,0]);
ub = double([Amwub,Aiwub,Aewub,r2smwub,r2siwub,r2sewub,fmwub,fiwub,100,100]);

if numMagn==numel(te) % magnitude fitting
    
else    % other fittings

    totalField0 = db0(:).';  totalFieldlb   = db0(:).'-25;                 totalFieldub    = db0(:).'+25;
    pini0       = pini(:).'; pinilb         = ones(size(pini0))*(-2*pi);    piniub          = ones(size(pini0))*2*pi;
    
    % set initial guess and fitting boundaries
    x0 = double([x0,totalField0,pini0]);
    lb = double([lb,totalFieldlb,pinilb]);
    ub = double([ub,totalFieldub,piniub]);
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
phiCycle = RF_phase_cycle(npulse,50);
for kfa=1:length(fa)
T3D_all{kfa} = PrecomputeT(phiCycle,d2r(fa(kfa)*b1));
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,npulse,T3D_all,numMagn,isWeighted,weightMethod,isNormCost,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,npulse,T3D_all,numMagn,isWeighted,weightMethod,isNormCost,DEBUG)
% capture all fitting parameters
Amw=x(1);   Aiw=x(2);   Aew=x(3);
t2smw=1/x(4); t2siw=1/x(5); t2sew=1/x(6);
% t1s=x(7);   
t1s = 1/4;
t1l = 1/0.41;
fmw=x(7);   fiw=x(8);  
freq_ew = 0; 
counter = 9;
% kls=x(counter); counter = counter+1; 
kmi=x(counter); counter = counter+1; 
kme=x(counter); counter = counter+1; 

if numMagn==numel(te) % magnitude fitting
    % no initial phase 
    pini=0;
    totalfield = zeros(size(fa));
else    % other fittings
    totalfield = x(counter:counter+length(fa)-1);
    pini       = x(counter+length(fa):end);
end



% simulate signal based on parameter input
% sHat = mwi_model_2T13T2scc_epgx(fa,te,tr,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1s,t1l,fmw,fiw,freq_ew,totalfield,pini,b1,kls,npulse,T3D_all);
% ss = mwi_model_3T1_ssSPGR(fa,tr,Amw,Aiw,Aew,t1s,t1l,t1l,b1);
fa = fa*b1;
% [~, ss] = mwi_model_ssSPGR_2T1(fa,tr,Amw,Aiw+Aew,t1s,t1l,[],kls);
[~,ss] = mwi_model_ssSPGR_3T1(fa,tr,Aiw,Amw,Aew,t1l,t1s,t1l,...
                                      kmi,kme);
% ss(3,:) = Aew / (Aiw+Aew) * ss(2,:);
% ss(2,:) = Aiw / (Aiw+Aew) * ss(2,:);
sHat = zeros(size(s));
for kfa = 1:length(fa)
%     sHat(kfa,:) = mwi_model_3cc_dimwi(te,ss(1,kfa),ss(2,kfa),ss(3,kfa),t2smw,t2siw,t2sew,fmw,fiw,freq_ew,totalfield(kfa),pini);
    sHat(kfa,:) = mwi_model_3cc_dimwi(te,ss(2,kfa),ss(1,kfa),ss(3,kfa),t2smw,t2siw,t2sew,fmw,fiw,freq_ew,totalfield(kfa),pini);
end


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

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
if isNormCost
    err = err ./ norm(abs(s(:)));
    % err = err ./ mean(abs(s(:)));
end

% if DEBUG then plots current fitting result
% Debug module
if DEBUG
    Debug_display(s,sHat,err,te,x,numMagn);
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
% check maximum iterations allowed
try algoPara2.maxIter = algoPara.maxIter;               catch; algoPara2.maxIter = 500; end
% check function tolerance
try algoPara2.fcnTol = algoPara.fcnTol;                 catch; algoPara2.fcnTol = 1e-5; end
% check step tolerance
try algoPara2.stepTol = algoPara.stepTol;               catch; algoPara2.stepTol = 1e-5; end
% check weighted sum of cost function
try algoPara2.weightMethod = algoPara.weightMethod;     catch; algoPara2.weightMethod = 'norm'; end
% check weighted sum of cost function
try algoPara2.isWeighted = algoPara.isWeighted;         catch; algoPara2.isWeighted = true; end
% check # of phase-corrupted echoes
try algoPara2.numMagn = algoPara.numMagn;               catch; algoPara2.numMagn = numel(imgPara.te); end
% check # of phase-corrupted echoes
try algoPara2.isNormCost = algoPara.isNormCost;         catch; algoPara2.isNormCost = true; end
% check user bounds and initial guesses
try algoPara2.userDefine.x0 = algoPara.userDefine.x0;   catch; algoPara2.userDefine.x0 = [];end
try algoPara2.userDefine.lb = algoPara.userDefine.lb;   catch; algoPara2.userDefine.lb = [];end
try algoPara2.userDefine.ub = algoPara.userDefine.ub;   catch; algoPara2.userDefine.ub = [];end
% check # of phase-corrupted echoes
try algoPara2.isInvivo = algoPara.isInvivo;             catch; algoPara2.isInvivo = true;   end
try algoPara2.npulse   = algoPara.npulse;               catch; algoPara.npulse    = 200;    end

%%%%%%%%%% 2. check data integrity %%%%%%%%%%
% check if the number of echo times matches with the data
if length(imgPara.te) ~= size(imgPara.img,4)
    error('The length of TE does not match with the 4th dimension of the image.');
end
% check if the number of flip angles matches with the data's 5th dim
if length(imgPara.fa) ~= size(imgPara.img,5)
    error('The length of FA does not match with the last dimension of the image.');
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

%%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
disp('Fitting options:');
fprintf('#RF pulse = %i\n',algoPara2.npulse);
fprintf('Max. iterations = %i\n',algoPara2.maxIter);
fprintf('Function tolerance = %.2e\n',algoPara2.fcnTol);
fprintf('Step tolerance = %.2e\n',algoPara2.stepTol);
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

disp('Cost function options:');
if algoPara2.isWeighted
    disp('Cost function weighted by echo intensity: True');
else
    disp('Cost function weighted by echo intensity: False');
end
if algoPara2.isNormCost
    disp('Cost function is normalised by signal intensity: True');
else
    disp('Cost function is normalised by signal intensity: False');
end

% determine the number of estimates
% numEst = 10; % basic setting has 6 estimates
numEst = 9; % basic setting has 6 estimates
if algoPara2.numMagn~=numel(imgPara2.te)
    numEst = numEst + length(imgPara2.fa) + 1; % total field and inital phase
end
algoPara2.numEst = numEst;

end

%%
function Debug_display(s,sHat,err,te,x,numMagn)
    global DEBUG_resnormAll
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    
    if mod(length(DEBUG_resnormAll),length(x)) == 1
    
    figure(99);

    if numMagn==numel(te)
        subplot(2,2,1);plot(abs(permute(s,[2 1])),'k^-');hold on;ylim([min(abs(s(:)))*0.9,max(abs(s(:)))*1.1]); title('Magnitude');
        plot(abs(permute(sHat,[2 1])),'x-.');hold off;
        subplot(2,2,2); 
        plot((abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'ro-.'); 
        title('residual');
        ha = subplot(2,2,3); pos = get(ha,'Position'); un = get(ha,'Units'); delete(ha)
        uitable('Data',x(:),'Units',un,'Position',pos);
        subplot(2,2,4);
        plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm')
        text(0.5,0.5,sprintf('resnorm=%e',sum(err(:).^2)),'Units','normalized');
    else
        subplot(2,3,1);
        plot(abs(permute(s,[2 1])),'k^-');hold on;
        plot(abs(permute(sHat,[2 1])),'x-.');hold off;
        ylim([min(abs(s(:)))*0.9,max(abs(s(:)))*1.1]);
        title('Magn.');
        
        subplot(2,3,2);
        plot(angle(permute(s,[2 1])),'k^-');hold on;
        plot(angle(permute(sHat,[2 1])),'x-');hold off;
        ylim([min(angle(s(:)))*0.9 max(angle(s(:)))*1.1]);
        title('Phase');
        
        subplot(2,3,3);
        plot(real(permute(s,[2 1])),'k^-');hold on;
        plot(imag(permute(s,[2 1])),'ks-');
        plot(real(permute(sHat,[2 1])),'bx-.');
        plot(imag(permute(sHat,[2 1])),'b*-.');hold off;
        ylim([min([real(s(:));imag(s(:))]),max([real(s(:));imag(s(:))])*1.1]);
        title('Real, Imaginary');
        
%         subplot(2,3,4);
%         plot(te(:).',real(permute(sHat(:)-s(:),[2 1])),'rx-.');hold on;
%         plot(te(:).',imag(permute(sHat(:)-s(:),[2 1])),'r*-.');hold off;
%         title('Residual');
        
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
end

%% progress display
function [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,mask)
previous_progress_percentage = floor(numFittedVoxel*100/numMaskedVoxel);

% update number of non zeros element in the current mask
numFittedVoxel = numFittedVoxel + nnz(mask);

current_progress_percentage = floor(numFittedVoxel*100/numMaskedVoxel);

if previous_progress_percentage ~= current_progress_percentage
    % delete previous progress
    for ii=1:length(progress)-1; fprintf('\b'); end
    % display current progress
    progress=sprintf('%d %%%%', floor(numFittedVoxel*100/numMaskedVoxel));
    fprintf(progress);
end
end