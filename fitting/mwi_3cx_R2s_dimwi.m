%% fitRes = mwi_3cx_R2s_dimwi(algoPara,imgPara)
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
% DIMWI model
% ===========
% 	.DIMWI.isVic     : true - use DIMWI to account for intra-axonal volume fraction; false - no DIMWI
% 	.DIMWI.isR2sEW   : true - use DIMWI to account for extra-axonal R2*; false - no DIMWI
% 	.DIMWI.isFreqMW  : true - use DIMWI to account for myelin water frequency; false - no DIMWI
% 	.DIMWI.isFreqIW  : true - use DIMWI to account for intra-axonal frequency; false - no DIMWI
%
% Advanced Starting point strategy
% ================================
%	.advancedStarting : 'default' - estimated S0 and T2s IEW from data
%                       'robust' - estimated S0 from data only
%                       otherwise fixed starting points for 3T
% imgPara       : structure array contains all image data
% Image-based
% ===========
% 	.img        : 4D image data, time in 4th dim
% 	.mask       : signal mask
% Acquisition parameters
% ======================
% 	.te         : echo times
% Recommended 
% ===========
% 	.fieldmap   : background field (default: 0)
% 	.pini       : initial phase 
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
% Output setting
% ==============
%   .output_dir      : directory to store final results (default:
%                      '/current_directory/mwi_results/')
%   .output_filename : output filename in text string (default:
%                      'mwi_results.mat')
%   .identifier      : temporary file identifier, a 8-digit code
%
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
%   .Freq_MW    : MW frequency, in Hz (MCR)
%   .Freq_IW    : IW frequency, in Hz (MCR)
%   .Freq_BKG   : Background frequency map, in Hz (complex-valued fitting)
%   .pini       : initial phase, in radian (complex-valued fitting)
%
% Description: Myelin water mapping using stanard or DIMWI model
% Chan, K.-S., Marques, J.P., 2020. Multi-compartment relaxometry and 
% diffusion informed myelin water imaging - promises and challenges of 
% new gradient echo myelin water imaging methods. Neuroimage 221, 117159. 
% https://doi.org/10.1016/j.neuroimage.2020.117159
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 February 2018
% Date modified: 16 August 2018
% Date modified: 6 March 2019
% Date modified: 22 October 2019
% Date modified: 29 October 2019
% Date modified: 17 Nov 2020
%
function fitRes = mwi_3cx_R2s_dimwi(algoPara,imgPara)

%%%%%%%%%% create directory for temporary results %%%%%%%%%%
default_output_filename = 'mwi_R2s_results.mat';
temp_prefix             = 'temp_mwi_R2s_';
[output_dir, output_filename, temp_filename, identifier] = set_up_output(imgPara,default_output_filename,temp_prefix);

%%%%%%%%%% log command window display to a text file %%%%%%%%%%
logFilename = fullfile(output_dir, ['run_mwi_' identifier '.log']);
logFilename = check_unique_filename_seqeunce(logFilename);
diary(logFilename)

% main
try
    
fprintf('Output directory                : %s\n',output_dir);
fprintf('Intermediate results identifier : %s\n',identifier);

disp('=====================================');
disp('Myelin water imaing: 3-pool T2* model');
disp('=====================================');

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

%%%%%%%%%% capture all image parameters %%%%%%%%%%
te    = double(imgPara.te);
data  = double(imgPara.img);
mask  = double(imgPara.mask);
fm    = double(imgPara.fieldmap);
pini  = double(imgPara.pini);

% basic info regarding data size
dims    = size(data);   dims = dims(1:3);
nVoxel  = numel(mask);
nTE     = length(te);

%%%%%%%%%% DIMWI %%%%%%%%%%
% check if fibre orientation is needed
if DIMWI.isVic
    icvf = double(imgPara.icvf);
else
    icvf = zeros(dims);
end
if DIMWI.isFreqMW || DIMWI.isFreqIW || DIMWI.isR2sEW
    b0dir = double(imgPara.b0dir);
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
    ff      = ones(dims);
    theta   = zeros(dims);
end

% store fixed parameters
fixParam.b0     = double(imgPara.b0);
fixParam.rho_mw = double(imgPara.rho_mw);
fixParam.E      = double(imgPara.E);
fixParam.x_i    = double(imgPara.x_i);
fixParam.x_a    = double(imgPara.x_a);
DIMWI.fixParam  = fixParam;

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

% normalise GRE image here
if isNormData
    
    [scaleFactor, data] = mwi_image_normalisation(data, mask);
    
%     tmp = abs(data(:,:,:,1));
% %     scaleFactor = norm(tmp(mask>0)) / sqrt(length(find(mask>0)));
else
    scaleFactor = 1;
end
% data = data / scaleFactor;

%%%%%%%%%% check advanced starting point strategy %%%%%%%%%%
advancedStarting = algoPara.advancedStarting;
if strcmpi(advancedStarting,'default') || strcmpi(advancedStarting,'robust')
    disp('Pre-estimate starting points using predefined model')
    % 3T
    if fixParam.b0 > 2.5 && fixParam.b0 < 3.5
        t2s_pre = [10e-3,60e-3]; % [T2sMW, T2sIEW] in second
    end
    
    switch advancedStarting
        case 'default'
            [s00,mwf0,t2siew0] = superfast_mwi_2m_standard_self(data,te,t2s_pre(1));
        case 'robust'
            [s00,mwf0] = superfast_mwi_2m_standard(data,te,t2s_pre);
    end
    s00 = sum(s00,4); % total water
    % also masked out problematic voxels detected by superfast method
    mask = and(mask,s00 ~= 0);
end

fitRes.mask_fitted = mask;
    
% reshape data to get only the masked data
mask    = reshape(mask,nVoxel,1);
% find masked voxels
ind     = find(mask~=0);
data   	= reshape(data, nVoxel,nTE);            data    = data(ind,:);
fm      = reshape(fm,   nVoxel,1);              fm      = fm(ind);
pini    = reshape(pini, nVoxel,1);              pini    = pini(ind);
icvf    = reshape(icvf, nVoxel,1);              icvf	= icvf(ind);
theta   = reshape(theta,nVoxel,size(theta,4));  theta   = theta(ind,:);
ff      = reshape(ff,   nVoxel,size(ff,4));     ff      = ff(ind,:);
if exist('s00','var');      s00     = reshape(s00,nVoxel,1);        s00     = s00(ind); end
if exist('mwf0','var');     mwf0    = reshape(mwf0,nVoxel,1);     	mwf0    = mwf0(ind); end
if exist('t2siew0','var');	t2siew0 = reshape(t2siew0,nVoxel,1);    t2siew0 = t2siew0(ind); end

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
% avoid too many voxels per batch which waits long to update progress
while nElement > 2000 && numBatch <=5000
    numBatch = numBatch + 50;
    nElement = floor(numMaskedVoxel/numBatch);
end
% avoid too few voxels per batch which reduces parallel computing efficiency
while nElement < 20 && numBatch ~= 1
    numBatch = round(numBatch/2);
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
    data_obj(kbat).data	= data(startInd:endInd,:);
    data_obj(kbat).fm	= fm(startInd:endInd);
    data_obj(kbat).pini	= pini(startInd:endInd);
    data_obj(kbat).icvf	= icvf(startInd:endInd);
    data_obj(kbat).theta= theta(startInd:endInd,:);
    data_obj(kbat).ff  	= ff(startInd:endInd,:);
    if exist('s00','var');      data_obj(kbat).s00      = s00(startInd:endInd); end
    if exist('mwf0','var');     data_obj(kbat).mwf0     = mwf0(startInd:endInd); end
    if exist('t2siew0','var');	data_obj(kbat).t2siew0	= t2siew0(startInd:endInd); end
end

%%%%%%%%%% initiate progress display %%%%%%%%%%
fprintf('%i voxel(s) to be fitted...\n',numMaskedVoxel);
if exist(temp_filename,'file')
    % restore progress
    disp('Restoring previous progress...')
    load(temp_filename);
    isRestore = true;
else
    % new progress
    fbat = 0;
    lastBatchEndTime = nElement *0.1; %roughly 0.1s for 1 voxel
    isRestore = false;
end
progress_display(numBatch,fbat,lastBatchEndTime,isRestore);
isRestore = false;

%%%%%%%%%% fitting main %%%%%%%%%%
if isParallel
    %%%%%%%%%% parfor loop %%%%%%%%%%
    for kbat = 1:numBatch
        
        data    = data_obj(kbat).data;
%         fm      = data_obj(kbat).fm;
%         pini    = data_obj(kbat).pini;
        icvf    = data_obj(kbat).icvf;
        theta   = data_obj(kbat).theta;
        ff      = data_obj(kbat).ff;
        
        for k = 1:size(data,1)
            initialGuess(k).pini0   = data_obj(kbat).pini(k);   % initial phase
            initialGuess(k).db0     = data_obj(kbat).fm(k);     % fieldmap
            if isfield(data_obj,'s00');     initialGuess(k).s00     = data_obj(kbat).s00(k);  end
            if isfield(data_obj,'mwf0');    initialGuess(k).mwf0    = data_obj(kbat).mwf0(k); end
            if isfield(data_obj,'t2siew0'); initialGuess(k).t2siew0 = data_obj(kbat).t2siew0(k); end
        end

        % start timer
        tic;
        % create an empty array for fitting results 
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        parfor k = 1:size(data,1)
            % T2*w
            s       = data(k,:);
%             db0     = fm(k);
%             pini0   = pini(k);
            icvf0   = icvf(k);
            theta0  = squeeze(theta(k,:));  theta0  = theta0(:);
            ff0     = squeeze(ff(k,:));     ff0     = ff0(:);
            initGuess = initialGuess(k);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,icvf0,theta0,ff0,initGuess,DIMWI,fitAlgor,userDefine,isInvivo,options,numEst,DEBUG);
        end
        lastBatchEndTime = toc;
        
        % Finished batch number
        fbat = kbat;
        
        res_obj(kbat).estimates    = estimates;
        res_obj(kbat).resnorm      = resnorm;
        res_obj(kbat).iterations   = iter;
        res_obj(kbat).exitflag     = exitflag;
        
        % display progress
        progress_display(numBatch,fbat,lastBatchEndTime,isRestore);
        save(temp_filename,'res_obj','fbat','lastBatchEndTime')

    end
else
    %%%%%%%%%% ordinary for loop %%%%%%%%%%
    for kbat=1:numBatch
        
        data    = data_obj(kbat).data;
%         fm      = data_obj(kbat).fm;
%         pini    = data_obj(kbat).pini;
        icvf    = data_obj(kbat).icvf;
        theta   = data_obj(kbat).theta;
        ff      = data_obj(kbat).ff;
        
        for k = 1:size(data,1)
            initialGuess(k).pini0   = data_obj(kbat).pini(k);   % initial phase
            initialGuess(k).db0     = data_obj(kbat).fm(k);     % fieldmap
            if isfield(data_obj,'s00');     initialGuess(k).s00     = data_obj(kbat).s00(k);  end
            if isfield(data_obj,'mwf0');    initialGuess(k).mwf0    = data_obj(kbat).mwf0(k); end
            if isfield(data_obj,'t2siew0'); initialGuess(k).t2siew0 = data_obj(kbat).t2siew0(k); end
        end
        
        % start timer
        tic;
        %%%%%%%%%% create an empty array for fitting results %%%%%%%%%%
        estimates = zeros(size(data,1),numEst);
        resnorm   = zeros(size(data,1),1);
        iter      = zeros(size(data,1),1);
        exitflag  = zeros(size(data,1),1);
        for k = 1:size(data,1)
            % T2*w
            s       = data(k,:);
%             db0     = fm(k);
%             pini0   = pini(k);
            icvf0   = icvf(k);
            theta0  = squeeze(theta(k,:));  theta0  = theta0(:);
            ff0     = squeeze(ff(k,:));     ff0     = ff0(:);
            initGuess = initialGuess(k);

            [estimates(k,:),resnorm(k),exitflag(k),iter(k)] = ...
                FitModel(s,te,icvf0,theta0,ff0,initGuess,DIMWI,fitAlgor,userDefine,isInvivo,options,numEst,DEBUG);
        end
        lastBatchEndTime = toc;
        
        % Finished batch number
        fbat = kbat;
        
        res_obj(kbat).estimates    = estimates;
        res_obj(kbat).resnorm      = resnorm;
        res_obj(kbat).iterations   = iter;
        res_obj(kbat).exitflag     = exitflag;
        
        % display progress
        progress_display(numBatch,fbat,lastBatchEndTime,isRestore);
        save(temp_filename,'res_obj','fbat','lastBatchEndTime')

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
estimates           = reshape(estimates,[dims,numEst]);
resnorm             = reshape(resnorm,dims);
iterations          = reshape(iterations,dims);
exitflag            = reshape(exitflag,dims);

%%%%%%%%%% Saving result %%%%%%%%%%
fitRes.estimates    = estimates;
fitRes.resnorm      = resnorm;
fitRes.iterations   = iterations;
fitRes.exitflag     = exitflag;
fitRes.mask_fitted  = and(exitflag ~= 99,fitRes.mask_fitted>0);
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

save(output_filename,'fitRes','-v7.3');
delete(temp_filename)
disp('Done!');

diary off

catch ME
    
    % close log file
    disp('There was an error! Please check the command window/error message file for more information.');
    diary off
    
    % open a new text file for error message
    errorMessageFilename = fullfile(output_dir, ['run_mwi_' identifier '.error']);
    errorMessageFilename = check_unique_filename_seqeunce(errorMessageFilename);
    fid = fopen(errorMessageFilename,'w');
    fprintf(fid,'The identifier was:\n%s\n\n',ME.identifier);
    fprintf(fid,'The message was:\n\n');
    msgString = getReport(ME,'extended','hyperlinks','off');
    fprintf(fid,'%s',msgString);
    fclose(fid);
    
    % rethrow the error message to command window
    rethrow(ME);
end

end

%% Setup lsqnonlin and fit with the default model
function [x,res,exitflag,iterations] = FitModel(s,te,icvf0,theta0,ff0,initialGuess,DIMWI,fitAlgor,userDefine,isInvivo,options,numEst,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end

s0 = abs(s(1));

try
b0 = DIMWI.fixParam.b0;
db0     = initialGuess.db0;
pini0   = initialGuess.pini0;
if isfield(initialGuess,'mwf0');    mwf     = min(initialGuess.mwf0,0.3); end   % maximum starting MWF is 30%
if isfield(initialGuess,'s00');     s0      = max(initialGuess.s00,s0); end     % start with highest S0 possible
if isfield(initialGuess,'t2siew0');	r2siew0	= 1/initialGuess.t2siew0; end       

if ~exist('r2siew0','var') || isnan(r2siew0) || isinf(r2siew0)
    r2siew0 = 18.5; % 3T
end

%%%%%%%%%% Step 1: determine and set initial guesses  %%%%%%%%%%
if isInvivo
    % range for in vivo, for 3T
    
%     r2smw0 = 100;	r2smwlb = 40;	r2smwub = 300;
%     r2siw0 = 16;	r2siwlb = 6;	r2siwub = 40;
%     r2sew0 = r2s0;	r2sewlb = 6;	r2sewub = 50;
    mwi_setup_initial_guess_3T_invivo;
else
    % range for ex vivo
    r2smw0 = 150;	r2smwlb = 40;	r2smwub = 300;
    r2siw0 = 20;	r2siwlb = 6;	r2siwub = 40;
    r2sew0 = r2s0;	r2sewlb = 6;	r2sewub = 50;
    if ~exist('mwf','var'); mwf = 0.15; end     % if not initial guess provided then set 15%
    v_ic = 0.8;
end

iwf = (1-mwf)*v_ic;

if fitAlgor.numMagn ~= numel(te)
    w_mb0 = db0*2*pi;	w_mblb = w_mb0-(15*b0)*2*pi;	w_mbub  = w_mb0+(15*b0)*2*pi;
    w_ib0 = db0*2*pi;	w_iblb = w_ib0-(8*b0)*2*pi;     w_ibub  = w_ib0+(3*b0)*2*pi;
else
    w_mb0 = 5*2*pi;	w_mblb = w_mb0-(15*b0)*2*pi;	w_mbub  = w_mb0+(15*b0)*2*pi;
    w_ib0 = 0*2*pi;	w_iblb = w_ib0-(8*b0)*2*pi;     w_ibub  = w_ib0+(3*b0)*2*pi;
end

% volume fraction of intra-axonal water
if DIMWI.isVic
    Amw0 = mwf*abs(s0);         Amylb = 0;	Amyub = abs(s0);
    Aie0 = (1-mwf)*abs(s0);     Aielb = 0;	Aieub = 1.5*abs(s0);
    
    x0 = double([Amw0 ,Aie0 ,r2smw0 ,r2siw0 ]);
    lb = double([Amylb,Aielb,r2smwlb,r2siwlb]);
    ub = double([Amyub,Aieub,r2smwub,r2siwub]);
else
    Amw0 = mwf*abs(s0);           Amylb = 0;	Amyub = abs(s0);
    Aiw0 = iwf*abs(s0);           Aiwlb = 0;	Aiwub = 1.5*abs(s0);
    Aew0 = (1-mwf-iwf)*abs(s0);  	Aewlb = 0;	Aewub = 1.5*abs(s0);
    
    x0 = double([Amw0 ,Aiw0 ,Aew0, r2smw0 ,r2siw0 ]);
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

catch 
    x           = zeros(1,numEst);
    res         = 0 ;
    exitflag    = 99;
    iterations  = 0;
end

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
    w = ones(size(s));
end
err = computeFiter(s,sHat,fitAlgor.numMagn,w);

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
% check maximum iterations allowed
try algoPara2.numBatch          = algoPara.numBatch;     	catch; algoPara2.numBatch = 50; end
% check normalised data before fitting
try algoPara2.isNormData        = algoPara.isNormData;  	catch; algoPara2.isNormData = true; end
% check maximum iterations allowed
try algoPara2.maxIter           = algoPara.maxIter;     	catch; algoPara2.maxIter = 500; end
% check function tolerance
try algoPara2.fcnTol            = algoPara.fcnTol;      	catch; algoPara2.fcnTol = 1e-5; end
% check step tolerance
try algoPara2.stepTol           = algoPara.stepTol;     	catch; algoPara2.stepTol = 1e-5; end
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
% check advanced starting points strategy
try algoPara2.advancedStarting	= algoPara.advancedStarting;catch; algoPara2.advancedStarting = 'default';end

%%%%%%%%%% 2. check data integrity %%%%%%%%%%
disp('-----------------------------');
disp('Checking input data integrity');
disp('-----------------------------');
% check if the number of echo times matches with the data
if length(imgPara.te) ~= size(imgPara.img,4)
    error('The length of TE does not match with the 4th dimension of the image.');
end
% check signal mask
disp('The following data is provided');
try
    imgPara2.mask = imgPara.mask;
    disp('Mask                          : True');
catch
    imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
    disp('Mask                          : false');
end
% check field map
try
    imgPara2.fieldmap = imgPara.fieldmap;
    disp('Field map                     : True');
catch
    if algoPara2.numMagn ~= length(imgPara2.te)
        imgPara2.fieldmap = angle(imgPara2.img(:,:,:,2)./imgPara2.img(:,:,:,1))/(2*pi*diff(imgPara2.te(1:2)));
    else
        imgPara2.fieldmap = zeros(size(imgPara2.mask));
    end
    disp('Field map                     : False');
end
% check initial phase map
try
    imgPara2.pini = imgPara.pini;
    disp('Initial phase input           : True');
catch
    imgPara2.pini = angle(exp(1i*(-2*pi*imgPara2.fieldmap*imgPara2.te(1)-angle(imgPara2.img(:,:,:,1)))));
    disp('Initial phase input           : False');
end
% check volume fraction of intra-axonal water
try
    imgPara2.icvf = imgPara.icvf;
    disp('Volume fraction of intra-axonal water input   : True');
catch
    imgPara2.icvf = zeros(size(imgPara2.mask));
    algoPara2.DIMWI.isVic = false;
    disp('Volume fraction of intra-axonal water input   : false');
end
% check fibre orientation map
try
    imgPara2.fo = imgPara.fo;
    disp('Fibre orientation input                       : True');
catch
    try
        imgPara2.theta = imgPara.theta;
        disp('Fibre orientation input                       : False');
        disp('Theta input                       : True');
    catch
        if algoPara2.DIMWI.isFreqMW || algoPara2.DIMWI.isFreqIW || algoPara2.DIMWI.isR2sEW
            error('Fibre orienation map or theta map is required for DIMWI');
        end
    end
end
% B0 direction
try
    imgPara2.b0dir = imgPara.b0dir;
    disp('B0dir input                   : True');
catch
    imgPara2.b0dir = [0,0,1];
    disp('B0dir input                   : false');
end
% field strength
try imgPara2.b0     = imgPara.b0;       catch; imgPara2.b0 = 3; end
try imgPara2.rho_mw = imgPara.rho_mw;   catch; imgPara2.rho_mw = 0.36/0.86; end
try imgPara2.x_i    = imgPara.x_i;    	catch; imgPara2.x_i = -0.1; end
try imgPara2.x_a    = imgPara.x_a;     	catch; imgPara2.x_a = -0.1; end
try imgPara2.E      = imgPara.E;     	catch; imgPara2.E = 0.02; end
disp('Input data is valid.')

%%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
disp('---------------');
disp('Fitting options');
disp('---------------');
if algoPara2.isNormData
    disp('GRE data is normalised before fitting');
else
    disp('GRE data is not normalised before fitting');
end
fprintf('Max. iterations    : %i\n',algoPara2.maxIter);
fprintf('Function tolerance : %.2e\n',algoPara2.fcnTol);
fprintf('Step tolerance     : %.2e\n',algoPara2.stepTol);
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
    disp('Initial guesses for in-vivo study');
else
    disp('Initial guesses for ex-vivo study');
end
disp('---------------------');
disp('Cost function options');
disp('---------------------');
% type of fitting
if algoPara2.numMagn==0
    disp('Fitting with complex data');
elseif algoPara2.numMagn==numel(imgPara2.te)
    disp('Fitting with magnitude data');
else
    fprintf('Fitting with %i magnitude data and %i complex data\n',algoPara2.numMagn,numel(imgPara2.te)-algoPara2.numMagn);
end
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
disp('------------------------------------');
disp('Diffusion informed MWI model options');
disp('------------------------------------');
if algoPara2.DIMWI.isVic
    disp('Volume fraction of intra-aonxal water will be used');
else
    disp('Volume fraction of intra-aonxal water will not be used');
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
disp('---------------------')
disp('Parameter to be fixed')
disp('---------------------')
disp(['Field strength (T)                       : ' num2str(imgPara2.b0)]);
disp(['Relative myelin water density            : ' num2str(imgPara2.rho_mw)]);
disp(['Myelin isotropic susceptibility (ppm)    : ' num2str(imgPara2.x_i)]);
disp(['Myelin anisotropic susceptibility (ppm)  : ' num2str(imgPara2.x_a)]);
disp(['Exchange term (ppm)                      : ' num2str(imgPara2.E)]);

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
function progress_display(numBatch,fbat,lastBatchEndTime,isRestore)
% function [progress,numFittedVoxel] = progress_display(progress,numMaskedVoxel,numFittedVoxel,newlyFittedVoxel)
% function [progress] = progress_display(progress,numMaskedVoxel,kind)

previous_progress_percentage    = 100*(fbat-1)/numBatch;
current_progress_percentage     = 100*fbat/numBatch;
% estimatedTime = current time for 1 batch * remaninig batches
estimatedTime = seconds(lastBatchEndTime*(numBatch-fbat));

if fbat ~= 0 && ~isRestore
    previous_progress = sprintf('\nProgress: %.2f %%%% \nEstimated remainaing time (HH:MM): %s\n', previous_progress_percentage,datestr(estimatedTime, 'HH:MM'));
    % delete previous progress
    for ii=1:length(previous_progress)-1; fprintf('\b'); end
end

% display current progress
progress = sprintf('\nProgress: %.2f %%%% \nEstimated remainaing time (HH:MM): %s\n', current_progress_percentage,datestr(estimatedTime, 'HH:MM'));
fprintf(progress);

end