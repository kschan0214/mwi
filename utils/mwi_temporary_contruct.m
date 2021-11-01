%% fitRes = mwi_3cx_2R1R2s_dimwi(algoPara,imgPara)
%
% Input
% --------------
% algoPara      : structure array contains all fitting algorithm specific
%                 parameters
%   .isInvivo   : initial starting points for in vivo study
%   .isFastEPG  : speed up the fitting by running a non EPG fitting first
%                 and use the results as the initial guesses of EPG fitting 
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
% Advanced Starting point strategy
% ================================
%	.advancedStarting : 'default' - estimated S0 and T2s IEW from data
%                       'robust' - estimated S0 from data only
%                       otherwise fixed starting points for 3T
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
% Output setting
% ==============
%   .output_dir      : directory to store final results (default:
%                      '/current_directory/mwi_results/')
%   .output_filename : output filename in text string (default:
%                      'mwi_results.mat')
%   .identifier      : temporary file identifier, a 8-digit code (optional)
%   .autoSave        : true (default); false - won't save anything after
%                      fitting
% Advanced Starting point strategy
% ================================
%	.advancedStarting : 'default' - estimated S0 and T2s IEW from data
%                       'robust' - estimated S0 from data only
%                       otherwise fixed starting points for 3T
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
%   .kiewm      : exchange rate, in s^-1 (fitted)
%   .Freq_MW    : MW frequency, in Hz (MCR)
%   .Freq_IW    : IW frequency, in Hz (MCR)
%   .Freq_BKG   : Background frequency map, in Hz (complex-valued fitting)
%   .pini       : initial phase, in radian (complex-valued fitting)
%
% Description: Myelin water mapping using MCR or MCR-DIMWI model
% Chan, K.-S., Marques, J.P., 2020. Multi-compartment relaxometry and 
% diffusion informed myelin water imaging - promises and challenges of 
% new gradient echo myelin water imaging methods. Neuroimage 221, 117159. 
% https://doi.org/10.1016/j.neuroimage.2020.117159
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 28 February 2018
% Date modified: 17 Nov 2020
% Date modified: 10 Aug 2021
%
function fitRes = mwi_temporary_contruct(algoPara,imgPara,res_obj)

disp('========================================');
disp('Testing');
disp('========================================');


% main
%%%%%%%%%% validate algorithm and image parameters %%%%%%%%%%
[algoPara,imgPara] = CheckAndSetDefault(algoPara,imgPara);


%%%%%%%%%% capture all algorithm parameters %%%%%%%%%%
numBatch   = algoPara.numBatch;
numEst     = algoPara.numEst;


%%%%%%%%%% capture all image parameters %%%%%%%%%%
te    = double(imgPara.te);
tr    = double(imgPara.tr);
fa    = double(imgPara.fa);
b0    = double(imgPara.b0);
data  = double(imgPara.img);
mask  = double(imgPara.mask);
b1map = double(imgPara.b1map);

% EPG-X related parameters
epgx.npulse       = algoPara.npulse;
epgx.rfphase      = algoPara.rfphase;
epgx.isExchange   = algoPara.isExchange;
epgx.isT1mw       = algoPara.isT1mw;
epgx.isEPG        = algoPara.isEPG;
epgx.T1mw         = algoPara.T1mw;
epgx.rho_mw       = double(imgPara.rho_mw);

% basic info regarding data size
dims    = size(data);   dims = dims(1:3);
nVoxel  = numel(mask);

%%%%%%%%%% DIMWI %%%%%%%%%%
[DIMWI, icvf, ff, ~] = setup_DIMWI(imgPara, algoPara);


%%%%%%%%%% prepare input %%%%%%%%%%
% parfor could be slowed down by a large number of unmasked voxels
% compute a new mask based on input data
mask    = and(mask>0,squeeze(sum(ff,4))>0);
if DIMWI.isVic
    mask = and(mask>0,icvf>0);
end


%%%%%%%%%% check advanced starting point strategy %%%%%%%%%%
advancedStarting = algoPara.advancedStarting;
if strcmpi(advancedStarting,'default') || strcmpi(advancedStarting,'robust')
    disp('Pre-estimate starting points using predefined model')
    % 3T
    if b0> 2.5 && b0 < 3.5
        t2s_pre = [10e-3,60e-3];    % [T2sMW, T2sIEW] in second
        t1_pre  = [234e-3, 1];    	% [T1MW, IEW], in second
    end
    
    switch advancedStarting
        case 'default'
            [s00,~,~,~] = superfast_mwi_2m_mcr_self(data,te,fa,tr,t2s_pre,t1_pre(1),mask,b1map,'superfast');
        case 'robust'
            [s00,~] = superfast_mwi_2m_mcr(data,te,fa,tr,t2s_pre,t1_pre,mask,b1map);
    end
    s00 = sum(s00,4); % total water
    % also masked out problematic voxels detected by superfast method
    mask = and(mask,s00 ~= 0);
end
% Simple DESPOT1 T1 mapping
if ~exist('t1iew0','var')
    [~, m0] = despot1_mapping(permute(data(:,:,:,1,:),[1 2 3 5 4]),fa,tr,mask,b1map);
    if ~exist('s00','var')
        s00 = m0;
        % also masked out DESPOT1 problematic voxels
        mask = and(mask,s00 ~= 0);
    end
end


% % find masked voxels
ind     = find(mask~=0);

%%%%%%%%%% Batch processing preparation %%%%%%%%%%
% get final no. of batches and elements of each batch
[numBatch, nElement] = setup_batch_determine_size(mask,numBatch);

% set up data_obj for batch processing
data_obj = setup_batch_create_data_obj(data,    mask, numBatch, nElement, 'data');

%%%%%%%%%% parfor loop %%%%%%%%%%
for kbat = 1:numBatch

    data    = data_obj(kbat).data;
    
    % create an empty array for fitting results 
    estimates = zeros(size(data,1),numEst);

    try 
        res_obj(kbat).estimates    = res_obj(kbat).estimates;
    catch 
        res_obj(kbat).estimates    = estimates;
    end

end


% concat results from all batches
tmp1    = [];
for kbat=1:numBatch
    tmp1 = cat(1,tmp1,res_obj(kbat).estimates);
end

% reshape fitting results back to map
estimates           = zeros(nVoxel,numEst);
estimates(ind,:)    = tmp1;
estimates           = reshape(estimates,[dims,numEst]);

%%%%%%%%%% Saving result %%%%%%%%%%
fitRes.estimates    = estimates;
% fitRes.mask_fitted  = and(exitflag ~= 99,fitRes.mask_fitted>0);
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


disp('Done!');

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
% fast fitting when EPG enabled
try algoPara2.isFastEPG       	= algoPara.isFastEPG;    	catch; algoPara2.isFastEPG = false; end
% check number of batches for parfor 
try algoPara2.numBatch          = algoPara.numBatch;     	catch; algoPara2.numBatch = 100; end
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
% check advanced starting points strategy
try algoPara2.advancedStarting  = algoPara.advancedStarting;catch; algoPara2.advancedStarting = [];end

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
    imgPara2.fieldmap = zeros([size(imgPara2.mask) length(imgPara.fa)]);
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
try imgPara2.b0     = imgPara.b0;       catch; imgPara2.b0      = 3; end
try imgPara2.rho_mw = imgPara.rho_mw;   catch; imgPara2.rho_mw  = 0.43; end
try imgPara2.x_i    = imgPara.x_i;    	catch; imgPara2.x_i     = -0.1; end
try imgPara2.x_a    = imgPara.x_a;     	catch; imgPara2.x_a     = -0.1; end
try imgPara2.E      = imgPara.E;     	catch; imgPara2.E       = 0.02; end

try imgPara2.autoSave = imgPara.autoSave; catch; imgPara2.autoSave = true; end
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
if algoPara2.isFastEPG
    disp('Fast EPG approach is used');
    warning('It does not yield the exact result with default fitting with EPG');
else
    disp('Fast EPG approach is not used');
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
    disp('Exchange - to be fitted');
else
    disp('Exchange - no exchange');
end
if algoPara2.isT1mw
    disp('T1mw - to be fitted');
else
    disp('T1mw - fixed');
    fprintf('T1mw will be fixed as %.3f s\n',algoPara2.T1mw);
end
if algoPara2.isEPG
    disp('EPG - True');
    fprintf('No. of RF pulses to equilibrium: %i\n', algoPara2.npulse);
    fprintf('Initial RF phase               : %i\n', algoPara2.rfphase);
else
    disp('EPG - no EPG simulation');
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
