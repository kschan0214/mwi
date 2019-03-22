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
% imgPara.pini     : initial phase 
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
% Date modified: 16 August 2018
% Date modified: 6 March 2019
%
%
function fitRes = mwi_2cx_T2s(algoPara,imgPara)
disp('Myelin water imaing: Multi-echo-T2* model');

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara]=CheckAndSetDefault(algoPara,imgPara);

% get debug mode
DEBUG   = algoPara.DEBUG;

% capture all algorithm parameters
numMagn    = algoPara.numMagn;
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
stepTol    = algoPara.stepTol;
isWeighted = algoPara.isWeighted;
isParallel = algoPara.isParallel;
userDefine = algoPara.userDefine;
isInvivo   = algoPara.isInvivo;
numEst     = algoPara.numEst;

te    = imgPara.te;
data  = imgPara.img;
mask  = imgPara.mask;
fm    = imgPara.fieldmap;
pini  = imgPara.pini;

[ny,nx,nz,~,~] = size(data);

% set fitting options
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*numEst,...
    'StepTolerance',stepTol,'FunctionTolerance',fcnTol);

if DEBUG
    % if DEBUG is on then disables parallel computing
    isParallel = false;
else
    options.Display = 'off';
end

% create empty array for fitting results
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
                    % T2*w
                    s       = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0     = fm(ky,kx,kz);
                    pini0   = pini(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,pini0,numMagn,isWeighted,userDefine,isInvivo,options,DEBUG);
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
                    % T2*w
                    s       = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0     = fm(ky,kx,kz);
                    pini0   = pini(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,pini0,numMagn,isWeighted,userDefine,isInvivo,options,DEBUG);
                end
            end
            % display progress
            [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,mask(ky,:,kz));
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

b = [ones(length(te),1), -te(:)]\log(abs(s(:)));
r2s = b(2);
s0 = exp(b(1));
if s0<0
    s0=0;
end

% set initial guesses
if isInvivo
    % in vivo reference
    Amw0   = 0.1*abs(s(1));     Amwlb   = 0;        Amwub   = 2*abs(s0);
    Aow0   = 0.9*abs(s(1));  	Aowlb   = 0;        Aowub   = 2*abs(s0);
    r2smw0 = 100;            	r2smwlb = 40;        r2smwub = 500;
%     t2siw0 = 64;                t2siwlb = 25;       t2siwub = 200;
    r2sow0 = r2s;                r2sowlb = 5;       r2sowub = 40;
else
    % ex vivo reference
    Amw0   = 0.15*abs(s(1));    Amwlb   = 0;        Amwub   = 2*abs(s0);
    Aow0   = 0.85*abs(s(1));  	Aowlb   = 0;        Aowub   = 2*abs(s0);
    r2smw0 = 10;            	r2smwlb = 1;        r2smwub = 25;
    r2sow0 = r2s;                r2sowlb = 25;       r2sowub = 200;
end
    
% set initial guess and fitting boundaries
x0 = double([Amw0 ,Aow0 ,r2smw0 ,r2sow0 ]);
lb = double([Amwlb,Aowlb,r2smwlb,r2sowlb]);
ub = double([Amwub,Aowub,r2smwub,r2sowub]);
if numMagn==numel(te) % magnitude fitting
    fmwbg0   = 5;           fmwbglb   = -75;    fmwbgub   = 75;
%     fowbg0   = 0;          	fiwbglb   = -25;    fiwbgub   = 25;

    % extra parameters
%     x0 = double([x0 ,fmwbg0 ,fowbg0]);
%     lb = double([lb,fmwbglb,fiwbglb]);
%     ub = double([ub,fmwbgub,fiwbgub]);
    x0 = double([x0 ,fmwbg0 ]);
    lb = double([lb,fmwbglb]);
    ub = double([ub,fmwbgub]);
else    % other fittings
    fmwbg0   = db0+5;         fmwbglb   = db0-75; 	fmwbgub   = db0+75;
    fowbg0   = db0;         fiwbglb   = db0-25;  	fiwbgub   = db0+25;
    if isnan(pini0)
        pini0  = angle(exp(1i*(-2*pi*db0*te(1)-angle(s(1)))));
    end
    pinilb = -2*pi;         piniub=2*pi;

    % extra parameters
    x0 = double([x0,fmwbg0,fowbg0,pini0]);
    lb = double([lb,fmwbglb,fiwbglb,pinilb]);
    ub = double([ub,fmwbgub,fiwbgub,piniub]);
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

if DEBUG
    x0
end

% run fitting algorithm here
[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,numMagn,isWeighted,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,numMagn,isWeighted,DEBUG)
% distribute fitting parameters
Amw=x(1);    Aow=x(2);
t2smw=1/x(3);  t2siw=1/x(4);
fmwbg=x(5);  
if numMagn==numel(te) % magnitude fitting
    pini=0; fiwbg=0; 
else    % other fittings
    pini=x(7); fiwbg=x(6); 
end

% simulate signal based on parameter input
sHat = mwi_model_2cc(te,Amw,Aow,t2smw,t2siw,fmwbg,fiwbg,pini);
% sHat = mwi_model_3cc_nam2015(te,Amw,Aow,Aew,t2smw,t2siw,t2sew,fmwbg,fiwbg,fbg,pini);

if size(sHat,1) ~= size(s,1)
    sHat = sHat.';
end

% compute fitting residual
if isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
%     w = sqrt(abs(s)/norm(abs(s(:))));
    w = sqrt(abs(s)/sum(abs(s(:))));
else
    % compute the cost without weights (=same weights)
    w = sqrt(ones(size(s))/numel(s));
end
err = computeFiter(s,sHat,numMagn,w);

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
err = err ./ norm(abs(s(:)));
% err = err ./ mean(abs(s(:)));

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
% check step tolerance
try algoPara2.stepTol = algoPara.stepTol;               catch; algoPara2.stepTol = 1e-6; end
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
disp('Fitting options:');
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

% determine the number of estimates
numEst = 5; % basic setting has 6 estimates
if algoPara2.numMagn~=numel(imgPara2.te)
    numEst = numEst + 2; % total field and inital phase
end
algoPara2.numEst = numEst;

end

%% Info display for debug mode
function Debug_display(s,sHat,err,te,x,numMagn)
    global DEBUG_resnormAll
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    figure(99);

    if numMagn==numel(te)
        subplot(2,2,1);plot(te(:).',abs(permute(s(:),[2 1])),'k^-');hold on;ylim([min(abs(s(:)))*0.9,max(abs(s(:)))*1.1]); title('Magnitude');
        plot(te(:).',abs(permute(sHat(:),[2 1])),'x-.');hold off;
        subplot(2,2,2); 
        plot(te(:).',(abs(permute(sHat(:),[2 1]))-abs(permute(s(:),[2 1]))),'ro-.'); 
        title('residual');
        ha = subplot(2,2,3); pos = get(ha,'Position'); un = get(ha,'Units'); delete(ha)
        uitable('Data',x(:),'Units',un,'Position',pos);
        subplot(2,2,4);
        plot(DEBUG_resnormAll);xlabel('# iterations');ylabel('resnorm')
        text(0.5,0.5,sprintf('resnorm=%e',sum(err(:).^2)),'Units','normalized');
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
        ylim([min([real(s(:));imag(s(:))]),max([real(s(:));imag(s(:))])*1.1]);
        title('Real, Imaginary');
        
        subplot(2,3,4);
        plot(te(:).',real(permute(sHat(:)-s(:),[2 1])),'rx-.');hold on;
        plot(te(:).',imag(permute(sHat(:)-s(:),[2 1])),'r*-.');hold off;
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

function [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,mask)
% number of non zeros element in the current mask
numFittedVoxel = numFittedVoxel + nnz(mask);
% delete previous progress
for ii=1:length(progress)-1; fprintf('\b'); end
% display current progress
progress=sprintf('%d %%%%', floor(numFittedVoxel*100/numMaskedVoxel));
fprintf(progress);

end