%% fitRes = mwi_3cx_R2s_dimwi(algoPara,imgPara)
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
% Date modified: 22 October 2019
% Date modified: 29 October 2019
%
function fitRes = mwi_3mm_R2s(algoPara,imgPara)
disp('Myelin water imaing: 3-pool T2* model');

%%%%%%%%%% validate algorithm and image parameters %%%%%%%%%%
[algoPara,imgPara] = CheckAndSetDefault(algoPara,imgPara);

%%%%%%%%%% get debug mode %%%%%%%%%%
DEBUG   = algoPara.DEBUG;

%%%%%%%%%% capture all algorithm parameters %%%%%%%%%%
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
stepTol    = algoPara.stepTol;
isParallel = algoPara.isParallel;
numBatch   = algoPara.numBatch;
numEst     = algoPara.numEst;
isInvivo   = algoPara.isInvivo;
userDefine = algoPara.userDefine;

fitAlgor.isWeighted = algoPara.isWeighted;
fitAlgor.isNormCost = algoPara.isNormCost;

%%%%%%%%%% capture all image parameters %%%%%%%%%%
te    = double(imgPara.te);
data  = double(imgPara.img);
mask  = double(imgPara.mask);

[ny,nx,nz,~,~] = size(data);

%%%%%%%%%% set fitting options %%%%%%%%%%
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*numEst,...
    'StepTolerance',stepTol,'FunctionTolerance',fcnTol);

if DEBUG
    % if DEBUG is on then disables parallel computing
    isParallel = false;
else
    options.Display = 'off';
end

%%%%%%%%%% prepare input %%%%%%%%%%
% parfor could be slowed down by a large number of unmasked voxels
% compute a new mask based on input data
    
mask    = reshape(mask,numel(mask),1);
% find masked voxels
ind     = find(mask~=0);
data   	= reshape(data, numel(mask),size(data,4));  data    = data(ind,:);

numMaskedVoxel = length(ind);
nElement = floor(numMaskedVoxel/numBatch);
for kbat = 1:numBatch
    startInd = (kbat-1)*nElement+1;
    if kbat < numBatch
        endInd      = kbat*nElement;
    else
        endInd      = length(ind);
    end
    data_obj(kbat).data	= data(startInd:endInd,:);
end

%%%%%%%%%% initiate progress display %%%%%%%%%%
numFittedVoxel = 0;
fprintf('%i voxel(s) to be fitted...\n',numMaskedVoxel);
progress='0 %%';
fprintf('Progress: ');
fprintf(progress);

%%%%%%%%%% fitting main %%%%%%%%%%
if isParallel
    %%%%%%%%%% parfor loop %%%%%%%%%%
    for kbat = 1:numBatch
        
        data    = data_obj(kbat).data;
        
        % create an empty array for fitting results 
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        parfor k = 1:size(data,1)
            % T2*w
            s       = data(k,:);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,fitAlgor,userDefine,isInvivo,options,DEBUG);
            
        end
        
        data_obj(kbat).estimates    = estimates;
        data_obj(kbat).resnorm      = resnorm;
        data_obj(kbat).iterations   = iter;
        data_obj(kbat).exitflag     = exitflag;
        newlyFittedVoxel = size(data,1);
        % display progress
        [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,newlyFittedVoxel);

    end
else
    %%%%%%%%%% ordinary for loop %%%%%%%%%%
    for kbat=1:numBatch
        
        data    = data_obj(kbat).data;
        
        %%%%%%%%%% create an empty array for fitting results %%%%%%%%%%
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        for k = 1:size(data,1)
            % T2*w
            s       = data(k,:);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,fitAlgor,userDefine,isInvivo,options,DEBUG);
            
        end
        data_obj(kbat).estimates = estimates;
        data_obj(kbat).resnorm   = resnorm;
        data_obj(kbat).iterations   = iter;
        data_obj(kbat).exitflag     = exitflag;
        newlyFittedVoxel = size(data,1);
        % display progress
        [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,newlyFittedVoxel);

    end
end
fprintf('\n');

% concat results from all batches
tmp1    = [];
tmp2    = [];
tmp3    = [];
tmp4    = [];
for kbat=1:numBatch
    tmp1 = cat(1,tmp1,data_obj(kbat).estimates);
    tmp2 = cat(1,tmp2,data_obj(kbat).resnorm);
    tmp3 = cat(1,tmp3,data_obj(kbat).iterations);
    tmp4 = cat(1,tmp4,data_obj(kbat).exitflag);
end

% reshape fitting results back to map
estimates           = zeros(ny*nx*nz,numEst);
resnorm             = zeros(ny*nx*nz,1);
iterations          = zeros(ny*nx*nz,1);
exitflag            = zeros(ny*nx*nz,1);
estimates(ind,:)    = tmp1;
resnorm(ind)        = tmp2;
iterations(ind)     = tmp3;
exitflag(ind)       = tmp4;
estimates           = reshape(estimates,ny,nx,nz,numEst);
resnorm             = reshape(resnorm,ny,nx,nz);
iterations          = reshape(iterations,ny,nx,nz);
exitflag            = reshape(exitflag,ny,nx,nz);

%%%%%%%%%% Saving result %%%%%%%%%%
fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;
fitRes.iterations = iterations;
fitRes.exitflag  = exitflag;
fitRes.mask      = mask;

fitRes.S0_MW = estimates(:,:,:,1); 
fitRes.S0_IW = estimates(:,:,:,2); 
fitRes.S0_EW = estimates(:,:,:,3); 
fitRes.R2s_MW = estimates(:,:,:,4); 
fitRes.R2s_IW = estimates(:,:,:,5); 
fitRes.R2s_EW = estimates(:,:,:,6); 

end

%% Setup lsqnonlin and fit with the default model
function [x,res,exitflag,iterations] = FitModel(s,te,fitAlgor,userDefine,isInvivo,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end

%%%%%%%%%% Step 1: determine and set initial guesses  %%%%%%%%%%
if isInvivo
    % range for in vivo
    r2smw0 = 100;	r2smwlb = 40;	r2smwub = 300;
    r2siw0 = 16;	r2siwlb = 6;	r2siwub = 40;
    r2sew0 = 21;	r2sewlb = 6;	r2sewub = 40;
    mwf = 0.1;
    iwf = 0.6;
else
    % range for ex vivo
    r2smw0 = 150;	r2smwlb = 40;	r2smwub = 300;
    r2siw0 = 20;	r2siwlb = 6;	r2siwub = 40;
    r2sew0 = 30;	r2sewlb = 6;	r2sewub = 40;
    mwf = 0.15;
    iwf = 0.6;
end

% volume fraction of intra-axonal water
Amy0 = mwf*abs(s(1));           Amylb = 0;	Amyub = 2*abs(s(1));
Aiw0 = iwf*abs(s(1));           Aiwlb = 0;	Aiwub = 2*abs(s(1));
Aew0 = (1-mwf-iwf)*abs(s(1));  	Aewlb = 0;	Aewub = 2*abs(s(1));

x0 = double([Amy0 ,Aiw0 ,Aew0, r2smw0 ,r2siw0 ,r2sew0]);
lb = double([Amylb,Aiwlb,Aewlb,r2smwlb,r2siwlb,r2sewlb]);
ub = double([Amyub,Aiwub,Aewub,r2smwub,r2siwub,r2sewub]);


%%%%%%%%%% Step 2: in case of user defined guesses  %%%%%%%%%%
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

%%%%%%%%%% Step 3: run fitting algorithm  %%%%%%%%%%
[x,res,~,exitflag,output] = lsqnonlin(@(y)CostFunc(y,s,te,fitAlgor,DEBUG),x0,lb,ub,options);

iterations = output.iterations;
end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,fitAlgor,DEBUG)
%%%%%%%%%% capture fitting parameters %%%%%%%%%%
% signal intensity
Amw=x(1);   Aiw=x(2);	Aew = x(3); 

% T2*s
t2smw = 1/x(4); t2siw = 1/x(5); t2sew = 1/x(6);

freq_mw = 0;	freq_iw = 0;   freq_ew = 0;
fbg = 0;    
pini = 0;


%%%%%%%%%% simulate signal based on parameter input %%%%%%%%%%
sHat = mwi_model_3cc_dimwi(te,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini);

% make sure the simulated data has the same size as the input
if size(sHat,1) ~= size(s,1)
    sHat = sHat.';
end

%%%%%%%%%% compute fitting residual %%%%%%%%%%
if fitAlgor.isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
    w = sqrt(abs(s));
else
    % compute the cost without weights (=same weights)
%     w = sqrt(ones(size(s))/numel(s));
    w = ones(size(s));
end
err = computeFiter(s,sHat,length(te),w);

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
if fitAlgor.isNormCost
    err = err ./ norm(abs(s(:)));
    % err = err ./ mean(abs(s(:)));
end

% Debug module
if DEBUG
    Debug_display(s,sHat,err,te,x,length(te));
end

end

%% check and set default
function [algoPara2,imgPara2]=CheckAndSetDefault(algoPara,imgPara)

imgPara2 = imgPara;
algoPara2 = algoPara;

%%%%%%%%%% 1. check algorithm parameters %%%%%%%%%%
% check debug
try algoPara2.DEBUG             = algoPara.DEBUG;         	catch; algoPara2.DEBUG = false; end
% check parallel computing 
try algoPara2.isParallel        = algoPara.isParallel;   	catch; algoPara2.isParallel = false; end
% check maximum iterations allowed
try algoPara2.numBatch          = algoPara.numBatch;     	catch; algoPara2.numBatch = 50; end
% check maximum iterations allowed
try algoPara2.maxIter           = algoPara.maxIter;     	catch; algoPara2.maxIter = 500; end
% check function tolerance
try algoPara2.fcnTol            = algoPara.fcnTol;      	catch; algoPara2.fcnTol = 1e-6; end
% check step tolerance
try algoPara2.stepTol           = algoPara.stepTol;     	catch; algoPara2.stepTol = 1e-6; end
% check weighted sum of cost function
try algoPara2.isWeighted        = algoPara.isWeighted;  	catch; algoPara2.isWeighted = true; end
% check # of phase-corrupted echoes
try algoPara2.isInvivo          = algoPara.isInvivo;       	catch; algoPara2.isInvivo = true; end
% check # of phase-corrupted echoes
try algoPara2.isNormCost        = algoPara.isNormCost;      catch; algoPara2.isNormCost = true; end
% check user bounds and initial guesses
try algoPara2.userDefine.x0     = algoPara.userDefine.x0;   catch; algoPara2.userDefine.x0 = [];end
try algoPara2.userDefine.lb     = algoPara.userDefine.lb;   catch; algoPara2.userDefine.lb = [];end
try algoPara2.userDefine.ub     = algoPara.userDefine.ub;   catch; algoPara2.userDefine.ub = [];end

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

%%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
disp('Fitting options:');
fprintf('Max. iterations    = %i\n',algoPara2.maxIter);
fprintf('Function tolerance = %.2e\n',algoPara2.fcnTol);
fprintf('Step tolerance     = %.2e\n',algoPara2.stepTol);
% type of fitting
disp('Fitting magnitude model with magnitude data');

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
numEst = 6; % basic setting has 4 estimates
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

%% progress display
function [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,newlyFittedVoxel)
% function [progress] = progress_display(progress,numMaskedVoxel,kind)

previous_progress_percentage    = floor(numFittedVoxel*100/numMaskedVoxel);

% update number of non zeros element in the current mask
numFittedVoxel = numFittedVoxel + newlyFittedVoxel;

current_progress_percentage = floor(numFittedVoxel*100/numMaskedVoxel);

if previous_progress_percentage ~= current_progress_percentage
    % delete previous progress
    for ii=1:length(progress)-1; fprintf('\b'); end
    % display current progress
    progress=sprintf('%d %%%%', floor(numFittedVoxel*100/numMaskedVoxel));
    fprintf(progress);
end
end