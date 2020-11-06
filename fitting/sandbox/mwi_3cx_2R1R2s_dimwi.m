%% fitRes = mwi_3cx_2R1R2s_dimwi(algoPara,imgPara)
%
% Input
% --------------
% algoPara      : structure array contains all fitting algorithm specific
%                 parameters
%   .isInvivo   : initial starting points for in vivo study
%   .isParallel : parallel computing using parfor
%   .userDefine : user defined starting set of estimated parameters
% Fitting algorithm option
% ========================
% 	.DEBUG      : debug mode (default false)
%   .maxIter    : maximum iteration allows (default 500)
% 	.fcnTol     : function tolerance (default: 1e-5)
% 	.stepTol	: step tolerance
% Residual option
% ===============
% 	.numMagn    : 0 - fitting with complex-valued data; length(TE) - fitting with magnitude data
% 	.isWeighted : cost function weighted by echo intensity (default: True)
% 	.weightMethod  : Weighted with respect to the first echo of each flip angle acqusition
% T1 model
% ========
% 	.isExchange : true - include BM equation to account for exchange effect;
%                 false - no magnetisation exchange
% 	.isEPG      : true - using EPG simulation; 
%                 false - standard steady-state equation
% 	.npulse     : no. of RF pulses to reach steady state in EPG
% 	.rfphase    : RF spoiling phase in EPG simulation, in degree
% 	.isT1mw     : true - fit myelin T1; false - fixed myelin T1; 
% 	.T1mw       : myelin T1 to be fixed, in second
% DIMWI model
% ===========
% 	.DIMWI.isVic     : true - use DIMWI to account for intra-axonal volume fraction; false - no DIMWI
% 	.DIMWI.isR2sEW   : true - use DIMWI to account for extra-axonal R2*; false - no DIMWI
% 	.DIMWI.isFreqMW  : true - use DIMWI to account for myelin water frequency; false - no DIMWI
% 	.DIMWI.isFreqIW  : true - use DIMWI to account for intra-axonal frequency; false - no DIMWI
%
% imgPara       : structure array contains all image data
% Image-based
% ===========
% 	.img        : 5D image data, time in 4th dim, flip angle in 5th dim
% 	.mask       : signal mask
% Acquisition parameters
% ======================
% 	.te         : echo times
%   .fa         : flip angles
% Recommended 
% ===========
% 	.fieldmap   : background field (default: 0)
% 	.pini       : initial phase 
%   .b1map      : B1 map
% DIMWI
% =====
% 	.theta      : angle between fibre and B0 directions, range [0, pi/2];
%   .ff         : fibre fraction, [row, col, slice, (fibre)]
%   .fo         : fibre direction maps, [row, col, slice, fibre vector, (fibre)]
% 	.icvf       : intra-axonal/cellular water volume fraction, range [0,1]
% Fixed parameters used for DIMWI or EPGX
% =======================================
% 	.b0     	: field strength, in tesla
%   .b0dir      : B0 field direction, [x,y,z]
% 	.rho_mw    	: relative myelin water density
% 	.E      	: exchange effect in signal phase, in ppm
% 	.x_i      	: myelin isotropic susceptibility, in ppm
% 	.x_a      	: myelin anisotropic susceptibility, in ppm
%
% Output
% --------------
% fitRes        : structure array contains all results and fitting report
% 	.estimates  : all fitting estimates 
% 	.resnorm    : L2 norm of fitting residual
%   .iterations : no. of fitting iterations used
%   .exitflag   : exitflag of lsqnonlin;
%   .S0_MW      : Proton density weighted MW signal
%   .S0_IW      : Proton density weighted IW signal (MCR)
%   .S0_EW      : Proton density weighted EW signal (MCR)
%   .S0_IEW     : Proton density weighted IEW signal (MCR-DIMWI)
%	.R2s_MW     : MW R2*, in s^-1
%	.R2s_IW     : IW R2*, in s^-1
%   .R2s_EW 	: IW R2*, in s^-1 (MCR)
%   .T1_MW      : MW T1, in second (fitted)
%   .T1_IEW     : IEW T1, in second
%   .kiewm      : exchange rate, in s^-1
%   .Freq_MW    : MW frequency, in Hz (MCR)
%   .Freq_IW    : IW frequency, in Hz (MCR)
%   .Freq_BKG   : Background frequency map, in Hz (complex-valued fitting)
%   .pini       : initial phase, in radian (complex-valued fitting)
%
% Description: Myelin water mapping using MCR or MCR-DIMWI model
% Chan, K.-S., Marques, J.P., 2020. Multi-compartment relaxometry and 
% diffusion informed myelin water imaging ? promises and challenges of 
% new gradient echo myelin water imaging methods. Neuroimage 221, 117159. 
% https://doi.org/10.1016/j.neuroimage.2020.117159
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 February 2018
% Date modified: 5 Nov 2020
%
function fitRes = mwi_3cx_2R1R2s_dimwi(algoPara,imgPara)
disp('Myelin water imaing: MCR/MCR-DIMWI model');

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

% data fitting method related parameters
fitAlgor.numMagn    = algoPara.numMagn;
fitAlgor.isWeighted = algoPara.isWeighted;
fitAlgor.weightMethod = algoPara.weightMethod;
% fitAlgor.isNormCost = algoPara.isNormCost;

% EPG-X related parameters
epgx.npulse       = algoPara.npulse;
epgx.rfphase      = algoPara.rfphase;
epgx.isExchange   = algoPara.isExchange;
epgx.isT1mw       = algoPara.isT1mw;
epgx.isEPG        = algoPara.isEPG;
epgx.T1mw         = algoPara.T1mw;

%%%%%%%%%% capture all image parameters %%%%%%%%%%
te    = double(imgPara.te);
tr    = double(imgPara.tr);
fa    = double(imgPara.fa);
b0dir = double(imgPara.b0dir);
data  = double(imgPara.img);
mask  = double(imgPara.mask);
fm    = double(imgPara.fieldmap);
pini  = double(imgPara.pini);
b1map = double(imgPara.b1map);

% basic info regarding data size
dims    = size(mask);
nVoxel  = numel(mask);
nTE     = length(te);
nFA     = lengt(fa);

%%%%%%%%%% DIMWI %%%%%%%%%%
% check if intra-axonal water volume fraction is needed
if DIMWI.isVic
    icvf = double(imgPara.icvf);
else
    icvf = zeros(dims);
end
% check if fibre orientation is needed
if DIMWI.isFreqMW || DIMWI.isFreqIW || DIMWI.isR2sEW
    ff    = double(imgPara.ff); % fibre fraction
    try 
        fo    = double(imgPara.fo); % fibre orientation w.r.t. B0
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
    ff      = ones(dims);
    theta   = zeros(dims);
end

% store fixed parameters
fixParam.b0     = double(imgPara.b0);
fixParam.rho_mw = double(imgPara.rho_mw);
fixParam.E      = double(imgPara.E);
fixParam.x_i    = double(imgPara.x_i);
fixParam.x_a    = double(imgPara.x_a);
DIMWI.fixParam = fixParam;
epgx.rho_mw = fixParam.rho_mw;

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
fitRes.mask      = mask;

if isNormData
    tmp = max(abs(data(:,:,:,1,:)),[],5);
    scaleFactor = norm(tmp(mask>0)) / sqrt(length(find(mask>0)));
else
    scaleFactor = 1;
end
data = data / scaleFactor;

mask    = reshape(mask,numel(mask),1);
% find masked voxels
ind     = find(mask~=0);
data   	= reshape(data, nVoxel,nTE,nFA);            data    = data(ind,:,:);
fm      = reshape(fm,   nVoxel,size(fm,4));      	fm      = fm(ind,:);
pini    = reshape(pini, nVoxel,size(pini,4));     	pini    = pini(ind,:);
icvf    = reshape(icvf, nVoxel,1);                	icvf	= icvf(ind);
theta   = reshape(theta,nVoxel,size(theta,4));     	theta   = theta(ind,:);
ff      = reshape(ff,   nVoxel,size(ff,4));       	ff      = ff(ind,:);
b1map 	= reshape(b1map,nVoxel,1);                	b1map  	= b1map(ind);

% repackaging for progress reporting and temporary storage
numMaskedVoxel = length(ind);
% check how many voxels per current batch
nElement = floor(numMaskedVoxel/numBatch);
% If the no. of voxels less than 200 or no. batches then do it in one go
if numMaskedVoxel < numBatch || numMaskedVoxel < 200
    numBatch = 1;
    % update nElement
    nElement = floor(numMaskedVoxel/numBatch);
end
% avoid too few voxels per batch which reduces parallel computing efficiency
while nElement < 20 && numBatch ~= 1
    numBatch = round(numBatch/2);
    nElement = floor(numMaskedVoxel/numBatch);
end
% avoid too many voxels per batch which waits long to update progress
while nElement > 5000 || numBatch <=1000
    numBatch = numBatch + 50;
    nElement = floor(numMaskedVoxel/numBatch);
end
% get final no. of voxles per batch
nElement = floor(numMaskedVoxel/numBatch);
for kbat = 1:numBatch
    startInd = (kbat-1)*nElement+1;
    if kbat < numBatch
        endInd      = kbat*nElement;
    else
        endInd      = length(ind);
    end
    data_obj(kbat).data	= data(startInd:endInd,:,:);
    data_obj(kbat).fm	= fm(startInd:endInd,:);
    data_obj(kbat).pini	= pini(startInd:endInd,:);
    data_obj(kbat).icvf	= icvf(startInd:endInd);
    data_obj(kbat).theta= theta(startInd:endInd,:);
    data_obj(kbat).ff  	= ff(startInd:endInd,:);
    data_obj(kbat).b1map= b1map(startInd:endInd);
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
        b1map   = data_obj(kbat).b1map;
        
        % create an empty array for fitting results 
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        parfor k = 1:size(data,1)
            % T2*w
            s       = squeeze(data(k,:,:));
            db0     = squeeze(fm(k,:));
            pini0   = squeeze(pini(k,:));
            icvf0   = icvf(k);
            theta0  = squeeze(theta(k,:));  theta0  = theta0(:);
            ff0     = squeeze(ff(k,:));     ff0     = ff0(:);
            b10     = b1map(k);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,tr,fa,b10,icvf0,theta0,ff0,db0,pini0,epgx,DIMWI,fitAlgor,userDefine,isInvivo,options,DEBUG);
        end
        
        res_obj(kbat).estimates    = estimates;
        res_obj(kbat).resnorm      = resnorm;
        res_obj(kbat).iterations   = iter;
        res_obj(kbat).exitflag     = exitflag;
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
        b1map   = data_obj(kbat).b1map;
        
        % create an empty array for fitting results 
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        for k = 1:size(data,1)
            % T2*w
            s       = squeeze(data(k,:,:));
            db0     = squeeze(fm(k,:));
            pini0   = squeeze(pini(k,:));
            icvf0   = icvf(k);
            theta0  = squeeze(theta(k,:));  theta0  = theta0(:);
            ff0     = squeeze(ff(k,:));     ff0     = ff0(:);
            b10     = b1map(k);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,tr,fa,b10,icvf0,theta0,ff0,db0,pini0,epgx,DIMWI,fitAlgor,userDefine,isInvivo,options,DEBUG);
        end
        
        res_obj(kbat).estimates    = estimates;
        res_obj(kbat).resnorm      = resnorm;
        res_obj(kbat).iterations   = iter;
        res_obj(kbat).exitflag     = exitflag;
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
    tmp1 = cat(1,tmp1,res_obj(kbat).estimates);
    tmp2 = cat(1,tmp2,res_obj(kbat).resnorm);
    tmp3 = cat(1,tmp3,res_obj(kbat).iterations);
    tmp4 = cat(1,tmp4,res_obj(kbat).exitflag);
end

% reshape fitting results back to map
estimates           = zeros(nVoxel,numEst);
resnorm             = zeros(nVoxel,1);
iterations          = zeros(nVoxel,1);
exitflag            = zeros(nVoxel,1);
estimates(ind,:)    = tmp1;
resnorm(ind)        = tmp2;
iterations(ind)     = tmp3;
exitflag(ind)       = tmp4;
estimates           = reshape(estimates,dim(1),dim(2),dim(3),numEst);
resnorm             = reshape(resnorm,dim(1),dim(2),dim(3));
iterations          = reshape(iterations,dim(1),dim(2),dim(3));
exitflag            = reshape(exitflag,dim(1),dim(2),dim(3));

%%%%%%%%%% Saving result %%%%%%%%%%
fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;
fitRes.iterations= iterations;
fitRes.exitflag  = exitflag;
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
if epgx.isT1mw
    fitRes.T1_MW = estimates(:,:,:,counter); counter= counter + 1;
end
fitRes.T1_IEW = estimates(:,:,:,counter); counter= counter + 1;
if epgx.isExchange
    fitRes.kiewm = estimates(:,:,:,counter); counter= counter + 1;
end
if ~DIMWI.isFreqMW
    fitRes.Freq_MW = estimates(:,:,:,counter)/(2*pi); counter= counter + 1;
end
if ~DIMWI.isFreqIW
    fitRes.Freq_IW = estimates(:,:,:,counter)/(2*pi); counter= counter + 1;
end
if fitAlgor.numMagn~=numel(te)
    fitRes.Freq_BKG = estimates(:,:,:,counter:counter+length(fa)-1)/(2*pi); 
    fitRes.pini = estimates(:,:,:,counter+length(fa):end);
end

end

%% Setup lsqnonlin and fit with the default model
function [x,res,exitflag,iterations] = FitModel(s,te,tr,fa,b10,icvf0,theta0,ff0,db0,pini0,epgx,DIMWI,fitAlgor,userDefine,isInvivo,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end
[t10,rho0] = DESPOT1(abs(s(1,:)),fa,tr,'b1',b10);

if t10<0 || t10 > 10 || isnan(t10) || isinf(t10)
    t10 = 1;
end
if rho0<0 || isinf(rho0) || isnan(rho0)
    rho0 = max(abs(s(:)))*10;
end

b0 = DIMWI.fixParam.b0;

%%%%%%%%%% Step 1: determine and set initial guesses  %%%%%%%%%%
if isInvivo
    % range for in vivo
    r2smw0 = 100;	r2smwlb = 40;	r2smwub = 300;
    r2siw0 = 16;	r2siwlb = 6;	r2siwub = 40;
    r2sew0 = 21;	r2sewlb = 6;	r2sewub = 50;
    mwf = 0.1;
    iwf = 0.6;
    t1s0   = 234e-3;  	t1slb   = 50e-3;  	t1sub   = 650e-3;
    t1l0   = t10;     	t1llb   = 500e-3; 	t1lub   = t10+1;
    kls0   = 0;       	klslb   = 0;      	klsub   = 20;       % exchange rate from long T1 to short T1
else
    % range for ex vivo
    r2smw0 = 150;	r2smwlb = 40;	r2smwub = 300;
    r2siw0 = 20;	r2siwlb = 6;	r2siwub = 40;
    r2sew0 = 30;	r2sewlb = 6;	r2sewub = 40;
    mwf = 0.15;
    iwf = 0.6;
    t1s0   = 225e-3;  	t1slb   = 50e-3;  	t1sub   = 650e-3;
    t1l0   = t10;     	t1llb   = 500e-3; 	t1lub   = t10+1;
    kls0    = 0;       	klslb    = 0;      	klsub    = 6;       % exchange rate from long T1 to short T1
end

% volume fraction of intra-axonal water
if DIMWI.isVic
    Amy0 = mwf*rho0;            Amylb = 0;	Amyub = 1.5*rho0;
    Aie0 = (1-mwf)*rho0;        Aielb = 0;	Aieub = 1.5*rho0;
    
    x0 = [Amy0 ,Aie0 ,r2smw0 ,r2siw0 ];
    lb = [Amylb,Aielb,r2smwlb,r2siwlb];
    ub = [Amyub,Aieub,r2smwub,r2siwub];
else
    Amy0 = mwf*rho0;           Amylb = 0;	Amyub = 1.5*rho0;
    Aiw0 = iwf*rho0;           Aiwlb = 0;	Aiwub = 1.5*rho0;
    Aew0 = (1-mwf-iwf)*rho0;  	Aewlb = 0;	Aewub = 1.5*rho0;
    
    x0 = [Amy0 ,Aiw0 ,Aew0, r2smw0 ,r2siw0 ];
    lb = [Amylb,Aiwlb,Aewlb,r2smwlb,r2siwlb];
    ub = [Amyub,Aiwub,Aewub,r2smwub,r2siwub];
end
% R2* of extracellular water
if ~DIMWI.isR2sEW
    x0 = [x0,r2sew0];
    lb = [lb,r2sewlb];
    ub = [ub,r2sewub];
end
% T1 and exchange
if epgx.isT1mw
    x0 = [x0,t1s0];
    lb = [lb,t1slb];
    ub = [ub,t1sub];
end
x0 = [x0,t1l0];
lb = [lb,t1llb];
ub = [ub,t1lub];
if epgx.isExchange
    x0 = [x0,kls0];
    lb = [lb,klslb];
    ub = [ub,klsub];
end

% Angular frequency of myelin water
if ~DIMWI.isFreqMW
    w_mw0 = 5*2*pi;	w_mwlb = (-5*b0)*2*pi;	w_mwub  = (5+10*b0)*2*pi;
    
    x0 = double([x0,w_mw0]);
    lb = double([lb,w_mwlb]);
    ub = double([ub,w_mwub]);
end
% Angular frequency of intra-axonal water
if ~DIMWI.isFreqIW
    w_iw0 = -2*2*pi;	w_iwlb = (-2-8*b0)*2*pi;	w_iwub  = (-2+8*b0)*2*pi;
    
    x0 = double([x0,w_iw0]);
    lb = double([lb,w_iwlb]);
    ub = double([ub,w_iwub]);
end

% extra parameters that depended on fitting method
if fitAlgor.numMagn~=numel(te) % non magnitude fittingfalse

    w_bg0 = db0(:).'*2*pi;	w_bglb = (db0(:).'-8*b0)*2*pi;      w_bgub  = (db0(:).'+8*b0)*2*pi;
    pini00 = pini0(:).';    pinilb = ones(size(pini0))*(-2*pi);	piniub  = ones(size(pini0))*2*pi;
%     w_bg0 = db0*2*pi;	w_bglb = (db0-25)*2*pi;	w_bgub  = (db0+25)*2*pi;
%     if isnan(pini0)
%         pini0  = angle(exp(1i*(-2*pi*db0*te(1)-angle(s(1)))));
%     end
%     pinilb = -2*pi;         piniub = 2*pi;

    % extra parameters
    x0 = double([x0,w_bg0, pini00]);
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

if epgx.isExchange && epgx.isEPG
    % precompute EPG-X's transition matrix here for speed
    phiCycle = RF_phase_cycle(epgx.npulse,epgx.rfphase);
    for kfa=1:length(fa)
        T3D_all{kfa} = PrecomputeT(phiCycle,d2r(fa(kfa)*b10));
    end
    epgx.T3D_all = T3D_all;
end

%%%%%%%%%% Step 3: run fitting algorithm  %%%%%%%%%%
[x,res,~,exitflag,output] = lsqnonlin(@(y)CostFunc(y,s,te,tr,fa,b10,epgx,DIMWI,fitAlgor,DEBUG),x0,lb,ub,options);

iterations = output.iterations;
end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,tr,fa,b1,epgx,DIMWI,fitAlgor,DEBUG)
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
% T1
if epgx.isT1mw
    t1s = x(counter);   counter = counter + 1;
else
    t1s = epgx.T1mw; 
end
t1l = x(counter);   counter = counter + 1;
if epgx.isExchange
    kls = x(counter);   counter = counter + 1;
else
    kls = 0;
end

% frequency shifts
if ~DIMWI.isFreqMW
    freq_mw = x(counter)/(2*pi); counter = counter + 1;
else
    freq_mw = 0;
end
if ~DIMWI.isFreqIW
    freq_iw = x(counter)/(2*pi); counter = counter + 1;
else
    freq_iw = 0;
end

% external scanner effects
if fitAlgor.numMagn==numel(te) % magnitude fitting
    fbg = zeros(size(fa));                          pini=0;
else    % other fittings
    fbg=x(counter:counter+length(fa)-1)/(2*pi);     pini=x(counter+length(fa):end);
end

freq_ew = 0;

%%%%%%%%%% simulate signal based on parameter input %%%%%%%%%%
sHat = zeros([size(s) length(DIMWI.theta)]);
DIMWI_curr = DIMWI;
for kfo = 1:length(DIMWI.theta)
    
    DIMWI_curr.theta = DIMWI.theta(kfo);
    
    sHat(:,:,kfo) = mwi_model_2T13T2scc_dimwi(te,tr,fa,b1,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1s,t1l,kls,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI_curr,epgx);
end
sHat = sum(bsxfun(@times,sHat,permute(DIMWI.ff(:),[2 3 1])),3);

%%%%%%%%%% compute fitting residual %%%%%%%%%%
% compute fitting residual
if fitAlgor.isWeighted
    switch fitAlgor.weightMethod
        case 'norm'
            % weights using echo intensity, as suggested in Nam's paper
            w = sqrt(abs(s));
%            w =  abs(s)/sum(abs(s(:)));
        case '1stEcho'
            % weights using the 1st echo intensity of each flip angle
            w = bsxfun(@rdivide,abs(s),abs(s(1,:)));
    end

else
    % compute the cost without weights
    w = ones(size(s));
end

err = computeFiter(s.',sHat.',fitAlgor.numMagn,w.');

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

imgPara2    = imgPara;
algoPara2   = algoPara;

%%%%%%%%%% 1. check algorithm parameters %%%%%%%%%%
% check debug
try algoPara2.DEBUG             = algoPara.DEBUG;         	catch; algoPara2.DEBUG = false; end
% check parallel computing 
try algoPara2.isParallel        = algoPara.isParallel;   	catch; algoPara2.isParallel = false; end
% check number of batches for parfor 
try algoPara2.numBatch          = algoPara.numBatch;     	catch; algoPara2.numBatch = 50; end
% check maximum iterations allowed
try algoPara2.maxIter           = algoPara.maxIter;     	catch; algoPara2.maxIter = 200; end
% check function tolerance
try algoPara2.fcnTol            = algoPara.fcnTol;      	catch; algoPara2.fcnTol = 1e-5; end
% check step tolerance
try algoPara2.stepTol           = algoPara.stepTol;     	catch; algoPara2.stepTol = 1e-5; end
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
% EPG-X
try algoPara2.isExchange        = algoPara.isExchange;      catch; algoPara2.isExchange = false; end
try algoPara2.isT1mw            = algoPara.isT1mw;          catch; algoPara2.isT1mw     = false; end
try algoPara2.T1mw              = algoPara.T1mw;            catch; algoPara2.T1mw       = 234e-3; end % Du et al. 2014 PLOS ONE
try algoPara2.isEPG             = algoPara.isEPG;           catch; algoPara2.isEPG      = false; end
try algoPara2.rfphase           = algoPara.rfphase;       	catch; algoPara2.rfphase    = 50; end
try algoPara2.npulse            = algoPara.npulse;       	catch; algoPara2.npulse     = 50; end

%%%%%%%%%% 2. check data integrity %%%%%%%%%%
disp('-----------------------');
disp('Checking data integrity');
disp('-----------------------');
% check if the number of echo times matches with the data
if length(imgPara.te) ~= size(imgPara.img,4)
    error('The length of TE does not match with the 4th dimension of the image.');
end
if length(imgPara.fa) ~= size(imgPara.img,5)
    error('The length of flip angle does not match with the 5th dimension of the image.');
end
% check signal mask
try
    imgPara2.mask = imgPara.mask;
    disp('Mask input                : True');
catch
    imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
    disp('Mask input                : false');
end
% check field map
try
    imgPara2.fieldmap = imgPara.fieldmap;
    disp('Field map input           : True');
catch
    imgPara2.fieldmap = zeros(size(imgPara2.mask));
    disp('Field map input           : False');
end
% check initial phase map
try
    imgPara2.pini = imgPara.pini;
    disp('Initial phase input       : True');
catch
    imgPara2.pini = ones(size(imgPara2.mask))*nan;
    disp('Initial phase input       : False');
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
disp('Input data is valid.')

%%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
disp('--------------');
disp('Fitting option');
disp('--------------');
if algoPara2.isNormData
    disp('GRE data is normalised before fitting');
else
    disp('GRE data is not normalised before fitting');
end
fprintf('Max. iterations    : %i\n',algoPara2.maxIter);
fprintf('Function tolerance : %.2e\n',algoPara2.fcnTol);
fprintf('Step tolerance     : %.2e\n',algoPara2.stepTol);
% type of fitting
if algoPara2.numMagn==0
    disp('Fitting with complex data');
elseif algoPara2.numMagn==numel(imgPara2.te)
    disp('Fitting with magnitude data');
else
    fprintf('Fitting with %i magnitude data and %i complex data\n',algoPara2.numMagn,numel(imgPara2.te)-algoPara2.numMagn);
end
% initial guess and fitting bounds
if isempty(algoPara2.userDefine.x0)
    disp('Default initial guess     : True');
else
    disp('Default initial guess     : False');
end
if isempty(algoPara2.userDefine.lb)
    disp('Default lower bound       : True');
else
    disp('Default lower bound       : False');
end
if isempty(algoPara2.userDefine.ub)
    disp('Default upper bound       : True');
else
    disp('Default upper bound       : False');
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
    disp(['Weighting method: ' algoPara2.weightMethod]);
else
    disp('Cost function weighted by echo intensity: False');
end
% if algoPara2.isNormCost
%     disp('Cost function is normalised by signal intensity: True');
% else
%     disp('Cost function is normalised by signal intensity: False');
% end

disp('--------------------------');
disp('Multi-compartment T1 model');
disp('--------------------------');
if algoPara2.isExchange
    disp('Exchange - True');
else
    disp('Exchange - False');
end
if algoPara2.isT1mw
    disp('T1mw - True');
else
    disp('T1mw - False');
    fprintf('T1mw will be fixed as %.3f s\n',algoPara2.T1mw);
end
if algoPara2.isEPG
    disp('EPG - True');
    fprintf('No. of RF pulses to equilibrium: %i\n', algoPara2.npulse);
    fprintf('Initial RF phase               : %i\n', algoPara2.rfphase);
else
    disp('EPG - False');
end

disp('------------------------------------');
disp('Diffusion informed MWI model options');
disp('------------------------------------');
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
    disp('R2*       - extra-cellular water estimated by HCFM');
else
    disp('R2*       - extra-cellular water to be fitted');
end

disp('-------------------------------')
disp('Parameter to be fixed for DIMWI')
disp('-------------------------------')
disp(['Field strength (T)                       : ' num2str(imgPara2.b0)]);
disp(['B0 direction(x,y,z)                      : ' num2str(imgPara2.b0dir(:)')]);
disp(['Relative myelin water density            : ' num2str(imgPara2.rho_mw)]);
disp(['Myelin isotropic susceptibility (ppm)    : ' num2str(imgPara2.x_i)]);
disp(['Myelin anisotropic susceptibility (ppm)  : ' num2str(imgPara2.x_a)]);
disp(['Exchange term (ppm)                      : ' num2str(imgPara2.E)]);

% determine the number of estimates
numEst = 5; % basic setting has 6 estimates
if algoPara2.isExchange
    numEst = numEst + 1;
end
if algoPara2.isT1mw
    numEst = numEst + 1;
end
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
    numEst = numEst + size(imgPara2.fieldmap,4) + size(imgPara2.pini,4); % total field and inital phase
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
            plot(te(:).',abs(permute(sHat,[2 1])),'x-.');hold off;
            ylim([min(abs(s(:)))*0.9,max(abs(s(:)))*1.1]);
            title('Magn.');

            subplot(2,3,2);
            plot(te(:).',angle(permute(s,[2 1])),'k^-');hold on;
            plot(te(:).',angle(permute(sHat,[2 1])),'x-');hold off;
            ylim([min(angle(s(:)))-abs(min(angle(s(:))))*0.1, max(angle(s(:)))+abs(max(angle(s(:))))*0.1]);
            title('Phase');

            subplot(2,3,3);
            plot(te(:).',real(permute(s,[2 1])),'k^-');hold on;
            plot(te(:).',imag(permute(s,[2 1])),'ks-');
            plot(te(:).',real(permute(sHat,[2 1])),'bx-.');
            plot(te(:).',imag(permute(sHat,[2 1])),'b*-.');hold off;
            ylim([min([real(s(:));imag(s(:))]),max([real(s(:));imag(s(:))])*1.1]);
            title('Real, Imaginary');

            subplot(2,3,4);
            plot(te(:).',real(permute(sHat-s,[2 1])),'rx-.');hold on;
            plot(te(:).',imag(permute(sHat-s,[2 1])),'r*-.');hold off;
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