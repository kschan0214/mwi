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
function fitRes = mwi_3cx_R2s_dimwi(algoPara,imgPara)
disp('Myelin water imaing: 3-pool T2* model');

%%%%%%%%%% validate algorithm and image parameters %%%%%%%%%%
[algoPara,imgPara] = CheckAndSetDefault(algoPara,imgPara);

%%%%%%%%%% get debug mode %%%%%%%%%%
DEBUG   = algoPara.DEBUG;

%%%%%%%%%% capture all algorithm parameters %%%%%%%%%%
isNormData = algoPara.isNormData;
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
stepTol    = algoPara.stepTol;
isParallel = algoPara.isParallel;
numBatch   = algoPara.numBatch;
numEst     = algoPara.numEst;
isInvivo   = algoPara.isInvivo;
DIMWI      = algoPara.DIMWI;
userDefine = algoPara.userDefine;

fitAlgor.numMagn    = algoPara.numMagn;
fitAlgor.isWeighted = algoPara.isWeighted;
% fitAlgor.isNormCost = algoPara.isNormCost;

%%%%%%%%%% capture all image parameters %%%%%%%%%%
te    = double(imgPara.te);
b0dir = double(imgPara.b0dir);
data  = double(imgPara.img);
mask  = double(imgPara.mask);
fm    = double(imgPara.fieldmap);
pini  = double(imgPara.pini);
% check if intra-axonal water volume fraction is needed
if DIMWI.isVic
    icvf = double(imgPara.icvf);
else
    icvf = zeros(size(mask));
end
% check if fibre orientation is needed
if DIMWI.isFreqMW || DIMWI.isFreqIW || DIMWI.isR2sEW
    ff    = double(imgPara.ff); % fibre fraction
    try 
        fo    = double(imgPara.fo); % fibre orientation
        theta = zeros(size(ff));
        for kfo = 1:size(fo,5)
            theta(:,:,:,kfo) = AngleBetweenV1MapAndB0(fo(:,:,:,:,kfo),b0dir);
        end
    catch 
        theta = double(imgPara.theta); % theta map
    end
    % normalise fibre fraction
    ff = bsxfun(@rdivide,ff,sum(ff,4));
    ff(isnan(ff)) = 0;
else
    ff      = ones(size(mask));
    theta   = zeros(size(mask));
end

% store fixed parameters
fixParam.b0     = double(imgPara.b0);
fixParam.rho_mw = double(imgPara.rho_mw);
fixParam.E      = double(imgPara.E);
fixParam.x_i    = double(imgPara.x_i);
fixParam.x_a    = double(imgPara.x_a);
DIMWI.fixParam = fixParam;

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
mask    = and(mask>0,squeeze(sum(ff,4))>0);
if DIMWI.isVic
    mask = and(mask>0,icvf>0);
end

if isNormData
    tmp = abs(data(:,:,:,1));
    scaleFactor = norm(tmp(mask>0)) / sqrt(length(find(mask>0)));
else
    scaleFactor = 1;
end
data = data / scaleFactor;
    
mask    = reshape(mask,numel(mask),1);
% find masked voxels
ind     = find(mask~=0);
data   	= reshape(data, numel(mask),size(data,4));  data    = data(ind,:);
fm      = reshape(fm,   numel(mask),1);             fm      = fm(ind);
pini    = reshape(pini, numel(mask),1);             pini    = pini(ind);
icvf    = reshape(icvf, numel(mask),1);             icvf	= icvf(ind);
theta   = reshape(theta,numel(mask),size(theta,4)); theta   = theta(ind,:);
ff      = reshape(ff,   numel(mask),size(ff,4));    ff      = ff(ind,:);

numMaskedVoxel = length(ind);
if numMaskedVoxel < numBatch
    numBatch = 1;
end
nElement = floor(numMaskedVoxel/numBatch);
for kbat = 1:numBatch
    startInd = (kbat-1)*nElement+1;
    if kbat < numBatch
        endInd      = kbat*nElement;
    else
        endInd      = length(ind);
    end
    data_obj(kbat).data	= data(startInd:endInd,:);
    data_obj(kbat).fm	= fm(startInd:endInd);
    data_obj(kbat).pini	= pini(startInd:endInd);
    data_obj(kbat).icvf	= icvf(startInd:endInd);
    data_obj(kbat).theta= theta(startInd:endInd,:);
    data_obj(kbat).ff  	= ff(startInd:endInd,:);
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
        fm      = data_obj(kbat).fm;
        pini    = data_obj(kbat).pini;
        icvf    = data_obj(kbat).icvf;
        theta   = data_obj(kbat).theta;
        ff      = data_obj(kbat).ff;
        
        % create an empty array for fitting results 
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        parfor k = 1:size(data,1)
            % T2*w
            s       = data(k,:);
            db0     = fm(k);
            pini0   = pini(k);
            icvf0   = icvf(k);
            theta0  = squeeze(theta(k,:));  theta0  = theta0(:);
            ff0     = squeeze(ff(k,:));     ff0     = ff0(:);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,icvf0,theta0,ff0,db0,pini0,DIMWI,fitAlgor,userDefine,isInvivo,options,DEBUG);
        end
        
        data_obj(kbat).estimates = estimates;
        data_obj(kbat).resnorm   = resnorm;
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
        fm      = data_obj(kbat).fm;
        pini    = data_obj(kbat).pini;
        icvf    = data_obj(kbat).icvf;
        theta   = data_obj(kbat).theta;
        ff      = data_obj(kbat).ff;
        
        %%%%%%%%%% create an empty array for fitting results %%%%%%%%%%
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        for k = 1:size(data,1)
            % T2*w
            s       = data(k,:);
            db0     = fm(k);
            pini0   = pini(k);
            icvf0   = icvf(k);
            theta0  = squeeze(theta(k,:));  theta0  = theta0(:);
            ff0     = squeeze(ff(k,:));     ff0     = ff0(:);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,icvf0,theta0,ff0,db0,pini0,DIMWI,fitAlgor,userDefine,isInvivo,options,DEBUG);
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
if DIMWI.isVic
    fitRes.estimates(:,:,:,1:2) = estimates(:,:,:,1:2)*scaleFactor;
else
    fitRes.estimates(:,:,:,1:3) = estimates(:,:,:,1:3)*scaleFactor;
end

counter = 1;
fitRes.S0_MW = estimates(:,:,:,counter)*scaleFactor; counter = counter+1;
if DIMWI.isVic
    fitRes.S0_IEW = estimates(:,:,:,counter)*scaleFactor; counter = counter+1;
else
    fitRes.S0_IW = estimates(:,:,:,counter)*scaleFactor; counter = counter+1;
    fitRes.S0_EW = estimates(:,:,:,counter)*scaleFactor; counter = counter+1;
end
fitRes.R2s_MW = estimates(:,:,:,counter); counter= counter + 1;
fitRes.R2s_IW = estimates(:,:,:,counter); counter= counter + 1;
if ~DIMWI.isR2sEW
    fitRes.R2s_EW = estimates(:,:,:,counter); counter= counter + 1;
end
if ~DIMWI.isFreqMW
    fitRes.Freq_MW = estimates(:,:,:,counter)/(2*pi); counter= counter + 1;
end
if ~DIMWI.isFreqIW
    fitRes.Freq_IW = estimates(:,:,:,counter)/(2*pi); counter= counter + 1;
end
if fitAlgor.numMagn~=numel(te)
    fitRes.Freq_BKG = estimates(:,:,:,counter)/(2*pi); counter= counter + 1;
    fitRes.pini = estimates(:,:,:,counter);
end
if ~DIMWI.isFreqMW && fitAlgor.numMagn~=numel(te)
    fitRes.Freq_MW = fitRes.Freq_MW - fitRes.Freq_BKG;
end
if ~DIMWI.isFreqIW && fitAlgor.numMagn~=numel(te)
    fitRes.Freq_IW = fitRes.Freq_IW - fitRes.Freq_BKG;
end

end

%% Setup lsqnonlin and fit with the default model
function [x,res,exitflag,iterations] = FitModel(s,te,icvf0,theta0,ff0,db0,pini0,DIMWI,fitAlgor,userDefine,isInvivo,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end

b0 = DIMWI.fixParam.b0;

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
if fitAlgor.numMagn ~= numel(te)
    w_mb0 = db0*2*pi;	w_mblb = (db0-25*b0)*2*pi;	w_mbub  = (db0+25*b0)*2*pi;
    w_ib0 = db0*2*pi;	w_iblb = (db0-8*b0)*2*pi;	w_ibub  = (db0+8*b0)*2*pi;
else
    w_mb0 = 5*2*pi;	w_mblb = (5-25*b0)*2*pi;	w_mbub  = (5+25*b0)*2*pi;
    w_ib0 = 0*2*pi;	w_iblb = (0-8*b0)*2*pi;     w_ibub  = (0+8*b0)*2*pi;
end

% volume fraction of intra-axonal water
if DIMWI.isVic
    Amy0 = mwf*abs(s(1));       Amylb = 0;	Amyub = 2*abs(s(1));
    Aie0 = (1-mwf)*abs(s(1));	Aielb = 0;	Aieub = 2*abs(s(1));
    
    x0 = double([Amy0 ,Aie0 ,r2smw0 ,r2siw0 ]);
    lb = double([Amylb,Aielb,r2smwlb,r2siwlb]);
    ub = double([Amyub,Aieub,r2smwub,r2siwub]);
else
    Amy0 = mwf*abs(s(1));           Amylb = 0;	Amyub = 2*abs(s(1));
    Aiw0 = iwf*abs(s(1));           Aiwlb = 0;	Aiwub = 2*abs(s(1));
    Aew0 = (1-mwf-iwf)*abs(s(1));  	Aewlb = 0;	Aewub = 2*abs(s(1));
    
    x0 = double([Amy0 ,Aiw0 ,Aew0, r2smw0 ,r2siw0 ]);
    lb = double([Amylb,Aiwlb,Aewlb,r2smwlb,r2siwlb]);
    ub = double([Amyub,Aiwub,Aewub,r2smwub,r2siwub]);
end
% R2* of extracellular water
if ~DIMWI.isR2sEW
    x0 = double([x0,r2sew0]);
    lb = double([lb,r2sewlb]);
    ub = double([ub,r2sewub]);
end

% Angular frequency of myelin water
if ~DIMWI.isFreqMW
    x0 = double([x0,w_mb0]);
    lb = double([lb,w_mblb]);
    ub = double([ub,w_mbub]);
end
% Angular frequency of intra-axonal water
if ~DIMWI.isFreqIW
    x0 = double([x0,w_ib0]);
    lb = double([lb,w_iblb]);
    ub = double([ub,w_ibub]);
end

% extra parameters that depended on fitting method
if fitAlgor.numMagn~=numel(te) % non magnitude fittingfalse

    w_bg0 = db0*2*pi;	w_bglb = (db0-8*b0)*2*pi;	w_bgub  = (db0+8*b0)*2*pi;
%     w_bg0 = db0*2*pi;	w_bglb = (db0-25)*2*pi;	w_bgub  = (db0+25)*2*pi;
    if isnan(pini0)
        pini0  = angle(exp(1i*(-2*pi*db0*te(1)-angle(s(1)))));
    end
    pinilb = -2*pi;         piniub = 2*pi;

    % extra parameters
    x0 = double([x0,w_bg0, pini0]);
    lb = double([lb,w_bglb,pinilb]);
    ub = double([ub,w_bgub,piniub]);
end

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

DIMWI.icvf  = icvf0;
DIMWI.theta = theta0;
DIMWI.ff    = ff0;

%%%%%%%%%% Step 3: run fitting algorithm  %%%%%%%%%%
[x,res,~,exitflag,output] = lsqnonlin(@(y)CostFunc(y,s,te,DIMWI,fitAlgor,DEBUG),x0,lb,ub,options);

iterations = output.iterations;
end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,DIMWI,fitAlgor,DEBUG)
%%%%%%%%%% capture fitting parameters %%%%%%%%%%
% signal intensity
Amw=x(1);    
if DIMWI.isVic
    Aie=x(2); 
    Aiw = Aie*DIMWI.icvf; Aew = Aie*(1-DIMWI.icvf);
    counter = 3;
else
    Aiw=x(2); Aew = x(3); 
    counter = 4;
end

% T2*s
t2smw=1/x(counter);         counter = counter + 1;
t2siw=1/x(counter);         counter = counter + 1;
if ~DIMWI.isR2sEW
    t2sew = 1/x(counter);   counter = counter + 1;
else
    t2sew = t2siw;
end

% frequency shifts
if ~DIMWI.isFreqMW
    freq_mwbg = x(counter)/(2*pi); counter = counter + 1;
else
    freq_mwbg = 0;
end
if ~DIMWI.isFreqIW
    freq_iwbg = x(counter)/(2*pi); counter = counter + 1;
else
    freq_iwbg = 0;
end

% external scanner effects
if fitAlgor.numMagn==numel(te) % magnitude fitting
    fbg=0;        pini=0;
else    % other fittings
    fbg=x(counter)/(2*pi);     pini=x(counter+1);
end

freq_mw = freq_mwbg - fbg;
freq_iw = freq_iwbg - fbg;
freq_ew = 0;

%%%%%%%%%% simulate signal based on parameter input %%%%%%%%%%
sHat = zeros([length(s) length(DIMWI.theta)]);
DIMWI_curr = DIMWI;
for kfo = 1:length(DIMWI.theta)
    
    DIMWI_curr.theta = DIMWI.theta(kfo);
    
    sHat(:,kfo) = mwi_model_3cc_dimwi(te,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI_curr);
    
end
sHat = sum(bsxfun(@times,sHat,DIMWI.ff(:).'),2);

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
err = computeFiter(s,sHat,fitAlgor.numMagn,w);

% % cost function is normalised with the norm of signal in order to provide
% % sort of consistence with fixed function tolerance
% if fitAlgor.isNormCost
%     err = err ./ norm(abs(s(:)));
%     % err = err ./ mean(abs(s(:)));
% end

% Debug module
if DEBUG
    Debug_display(s,sHat,err,te,x,fitAlgor.numMagn);
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
% check normalised data before fitting
try algoPara2.isNormData        = algoPara.isNormData;  	catch; algoPara2.isNormData = true; end
% check weighted sum of cost function
try algoPara2.isWeighted        = algoPara.isWeighted;  	catch; algoPara2.isWeighted = true; end
% check # of phase-corrupted echoes
try algoPara2.numMagn           = algoPara.numMagn;         catch; algoPara2.numMagn = numel(imgPara.te); end
% check # of phase-corrupted echoes
try algoPara2.isInvivo          = algoPara.isInvivo;       	catch; algoPara2.isInvivo = true; end
% % check # of phase-corrupted echoes
% try algoPara2.isNormCost        = algoPara.isNormCost;      catch; algoPara2.isNormCost = true; end
% check user bounds and initial guesses
try algoPara2.userDefine.x0     = algoPara.userDefine.x0;   catch; algoPara2.userDefine.x0 = [];end
try algoPara2.userDefine.lb     = algoPara.userDefine.lb;   catch; algoPara2.userDefine.lb = [];end
try algoPara2.userDefine.ub     = algoPara.userDefine.ub;   catch; algoPara2.userDefine.ub = [];end
% check hollow cylinder fibre model parameters
try algoPara2.DIMWI.isFreqMW	= algoPara.DIMWI.isFreqMW;	catch; algoPara2.DIMWI.isFreqMW = true;end
try algoPara2.DIMWI.isFreqIW	= algoPara.DIMWI.isFreqIW;	catch; algoPara2.DIMWI.isFreqIW = true;end
try algoPara2.DIMWI.isR2sEW     = algoPara.DIMWI.isR2sEW;	catch; algoPara2.DIMWI.isR2sEW  = false;end
try algoPara2.DIMWI.isVic       = algoPara.DIMWI.isVic;     catch; algoPara2.DIMWI.isVic    = false;end

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
% check initial phase map
try
    imgPara2.pini = imgPara.pini;
    disp('Initial phase input: True');
catch
    imgPara2.pini = ones(size(imgPara2.mask))*nan;
    disp('Initial phase input: False');
end
% check volume fraction of intra-axonal water
try
    imgPara2.icvf = imgPara.icvf;
    disp('Volume fraction of intra-axonal water input: True');
end
% check fibre orientation map
try
    imgPara2.fo = imgPara.fo;
    disp('Fibre orientation input: True');
catch
    try
        imgPara2.theta = imgPara.theta;
        disp('Fibre orientation input: False');
        disp('Theta input: True');
    catch
        if algoPara2.DIMWI.isFreqMW || algoPara2.DIMWI.isFreqIW || algoPara2.DIMWI.isR2sEW
            error('Fibre orienation map or theta map is required for DIMWI');
        end
    end
end
% B0 direction
try
    imgPara2.b0dir = imgPara.b0dir;
    disp('B0dir input: True');
catch
    imgPara2.b0dir = [0,0,1];
    disp('B0dir input: false');
end
% field strength
try imgPara2.b0     = imgPara.b0;       catch; imgPara2.b0 = 3; end
try imgPara2.rho_mw = imgPara.rho_mw;   catch; imgPara2.rho_mw = 0.43; end
try imgPara2.x_i    = imgPara.x_i;    	catch; imgPara2.x_i = -0.1; end
try imgPara2.x_a    = imgPara.x_a;     	catch; imgPara2.x_a = -0.1; end
try imgPara2.E      = imgPara.E;     	catch; imgPara2.E = 0.02; end
disp('Parameter to be fixed:')
disp('----------------------')
disp(['Field strength (T)                       : ' num2str(imgPara2.b0)]);
disp(['Relative myelin water density            : ' num2str(imgPara2.rho_mw)]);
disp(['Myelin isotropic susceptibility (ppm)    : ' num2str(imgPara2.x_i)]);
disp(['Myelin anisotropic susceptibility (ppm)  : ' num2str(imgPara2.x_a)]);
disp(['Exchange term (ppm)                      : ' num2str(imgPara2.E)]);


% % check intra-aonxal water volume fraction map
% try
%     imgPara2.v_ic = imgPara.v_ic;
%     algoPara2.DIMWI.isVic = true;
%     disp('Intra-axonal water volume fraction input: True');
% catch
%     imgPara2.v_ic = zeros(size(imgPara2.mask));
%     algoPara2.DIMWI.isVic = false;
%     disp('Intra-axonal water volume fraction input: False');
% end

%%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
disp('Fitting options:');
disp('----------------');
if algoPara2.isNormData
    disp('GRE data is normalised before fitting');
else
    disp('GRE data is not normalised before fitting');
end
fprintf('Max. iterations    = %i\n',algoPara2.maxIter);
fprintf('Function tolerance = %.2e\n',algoPara2.fcnTol);
fprintf('Step tolerance     = %.2e\n',algoPara2.stepTol);
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
% if algoPara2.isNormCost
%     disp('Cost function is normalised by signal intensity: True');
% else
%     disp('Cost function is normalised by signal intensity: False');
% end

disp('Diffusion informed MWI model options:');
if algoPara2.DIMWI.isVic
    disp('Volume fraction of intra-aonxal water is provided');
else
    disp('Volume fraction of intra-aonxal water is NOT provided');
end
if algoPara2.DIMWI.isFreqMW
    disp('Frequency - myelin water estimated by HCFM');
else
    disp('Frequency - myelin water to be fitted');
end
if algoPara2.DIMWI.isFreqIW
    disp('Frequency - intra-axonal water estimated by HCFM');
else
    disp('Frequency - intra-axonal water to be fitted');
end
if algoPara2.DIMWI.isR2sEW
    disp('R2* - extra-cellular water estimated by HCFM');
else
    disp('R2* - extra-cellular water to be fitted');
end

% determine the number of estimates
numEst = 4; % basic setting has 4 estimates
if ~algoPara2.DIMWI.isVic
    numEst = numEst + 1;
end
if ~algoPara2.DIMWI.isFreqMW
    numEst = numEst + 1;
end
if ~algoPara2.DIMWI.isFreqIW
    numEst = numEst + 1;
end
if ~algoPara2.DIMWI.isR2sEW
    numEst = numEst + 1;
end
if algoPara2.numMagn~=numel(imgPara2.te)
    numEst = numEst + 2; % total field and inital phase
end
algoPara2.numEst = numEst;

end

%% Info display for debug mode
function Debug_display(s,sHat,err,te,x,numMagn)
    global DEBUG_resnormAll
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    
    if mod(numel(DEBUG_resnormAll), numel(x)) == 1
    
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
            ylim([min(angle(s(:)))-abs(min(angle(s(:))))*0.1, max(angle(s(:)))+abs(max(angle(s(:))))*0.1]);
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