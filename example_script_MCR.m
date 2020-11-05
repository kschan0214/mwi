addpath(genpath('epgx'));   % https://github.com/mriphysics/EPG-X
addpath(genpath('despot1')) % https://github.com/kschan0214/despot1
addpath(genpath('mwi'));    % https://github.com/kschan0214/mwi
% depending on the Matlab version, additional tool for NIFTI I/O might be
% required. I am using
% https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% here but it can be changed to whatever works for you
clear;

%% loading data
gre_dir         = ''; % directory of MGRE magnitude and phase data
mask_dir        = ''; % directory of signal
b1_dir          = ''; % directory of B1 map
sepia_res_dir   = [gre_dir 'sepia/']; % I strongly recommand using SEPIA (or other means) to obtain the total field map for each flip angle acquisition prior MWI

mask_fn         = ''; % mask NFITI filename
b1_fn           = ''; % B1 NIFTI filename, registered to GRE space already

% acquisition parameters
TR = TR;            % repetition time, in second
TE = [TE1,TE2,TE3]; % all echo times, in second, order should be the same as the way image data ared loaded
fa = [FA1,FA2,FA3]; % all flip angle, in degree, order should be the same as the way image data ared loaded

% load mask and B1 map
mask_nii    = load_untouch_nii(fullfile(mask_dir, mask_fn));
mask        = double(mask_nii.img);
b1_nii      = load_untouch_nii(fullfile(b1_dir ,b1_fn));
b1          = double(b1_nii.img);

matrixSize = size(mask);

% img is a 5D complex-valued data, [row,col,slice,echo time,flip angle];
img         = zeros([matrixSize length(TE) length(fa)]);
totalField  = zeros([matrixSize length(fa)]);
% load variable flip angle MGRE data
for kfa = 1:length(fa)

    magn_fn         = ['']; % MGRE magnitude data filename of the n-th flip angle
    phase_fn        = ['']; % MGRE phase data filename of the n-th flip angle
    totalField_fn   = ['']; % MGRE total field filename of the n-th flip angle
%     real_fn         = ['']; % MGRE real data filename of the n-th flip angle
%     imaginary_fn    = ['']; % MGRE imaginary data filename of the n-th flip angle
    
    magn_nii        = load_untouch_nii(fullfile(gre_dir, magn_fn));
    phase_nii       = load_untouch_nii(fullfile(gre_dir, phase_fn));
    totalField_nii 	= load_untouch_nii(fullfile(sepia_res_dir, totalField_fn));
%     real_nii        = load_untouch_nii(fullfile(gre_dir, real_fn));
%     imaginary_nii   = load_untouch_nii(fullfile(gre_dir, imaginary_fn));
    
    % convert to complex-valued data
    img(:,:,:,:,kfa)   = magn_nii.img .* exp(1i*phase_nii.img);
%     img(:,:,:,:,kfa)   = real_nii.img + 1i*imaginary_nii.img;
    
    totalField(:,:,:,kfa) = totalField_nii.img; % total field, in Hz
end

% estimate inital phase casued by transmit/receive B1, in radian
pini = zeros([matrixSize length(fa)]);
for kfa = 1:length(fa)
    pini(:,:,:,kfa) = angle(img(:,:,:,1,kfa) ./ exp(1i*2*pi*totalField(:,:,:,kfa)*TE(1)));
end
pini = mean(pini,4);

%% setup parameters
%%%% Step 1: setup fitting algorithm parameters
algoParam.isInvivo      = true;     % using initial starting points for in vivo study
algoParam.isParallel    = false;    % no parallel computing, can be turn on when set to true
algoParam.DEBUG         = false;    % no DEBUG mode for development
% fitting option
algoParam.maxIter       = 200;      % maximum number of fitting iteration
algoParam.fcnTol        = 1e-4;     % stopping tolerance, optimal values often between 1e-4 and 1e-6
algoParam.stepTol       = 1e-4;     % stopping tolerance, optimal values often between 1e-4 and 1e-6
% residual option
algoParam.numMagn       = 0;        % 0: fitting with complex-valued data; length(TE): fitting with magnitude data
algoParam.isWeighted    = true;     % SNR weighting of the fitting cost 
algoParam.weightMethod  = '1stEcho';% Weighted with respect to the first echo of each flip angle acqusition
% T1 model
algoParam.isExchange    = true;     % include BM equation to account for exchange effect
algoParam.isEPG         = true;     % using EPG simulation, can be switch off when set to false
algoParam.npulse        = 50;       % no. of RF pulses to reach steady state
algoParam.rfphase       = 50;       % RF spoiling phase for EPG simulation, in degree
algoParam.isT1mw        = false;    % false: fixed myelin T1; true: estimate myelin T1 in fitting
algoParam.T1mw          = 234e-3;   % value of myelin T1 to be fixed, in second
% DIMWI model
algoParam.DIMWI.isVic       = false;    % not using DIMWI to account for intra-axonal volume fraction
algoParam.DIMWI.isR2sEW     = false;    % not using DIMWI to account for extra-axonal R2*
algoParam.DIMWI.isFreqMW    = false;    % not using DIMWI to account for myelin water frequency
algoParam.DIMWI.isFreqIW    = false;    % not using DIMWI to account for intra-axonal frequency
% fixed parameters
kappa_mw            = 0.36; % Jung, NI., myelin water density
kappa_iew           = 0.86; % Jung, NI., intra-/extra-axonal water density
imgParam.b0     	= 3;    % field strength, in tesla
imgParam.rho_mw    	= kappa_mw/kappa_iew; % relative myelin water density
imgParam.E      	= 0.02; % exchange effect in signal phase, in ppm
imgParam.x_i      	= -0.1; % myelin isotropic susceptibility, in ppm
imgParam.x_a      	= -0.1; % myelin anisotropic susceptibility, in ppm

%%%% Step 2: input image normalisation
% normalise image for stable fitting result
% recommandation: performs nomalisation outside the fitting function if you
% want to manually run the processing in multiple CPUs in parallel
% if you decide to normalise the data outside the function, you should
% apply the same scaling factor to the signal ampltiude results, see below
algoParam.isNormData    = false;     % no normalised data before data fitting
tmp         = max(abs(img(:,:,:,1,:)),[],5);
scaleFactor = norm(tmp(mask>0)) / sqrt(length(find(mask>0)));
% same operation can be performed by the following setting if you use the function built-in parallelisation: 
% algoParam.isNormData    = true;      % normalised data before data fitting

%%%% Step 3: Constructing input data for the function
imgParam.te         = TE;
imgParam.tr         = TR;
imgParam.fa         = fa;
imgParam.img        = img/scaleFactor;
imgParam.mask       = mask;
imgParam.fieldmap   = totalField;       % it is not neccessary to provide the total field for initial starting point, yet it can provide the fitting results
imgParam.pini       = pini;
imgParam.b1map      = b1;

%%%% Step 4: fitting core
tic
fitres = mwi_3cx_2R1R2s_dimwi(algoParam,imgParam);
toc

%%%% Step 5: Rescale data (optional)
% if you normalise the data outside the fitting function, then you should
% rescale your data now
% you can omit this step if you use the built-in normalisation
fitres.estimates(:,:,:,1:3) = fitres.estimates(:,:,:,1:3)*scaleFactor;
fitres.S0_MW                = fitres.S0_MW*scaleFactor;
fitres.S0_IW                = fitres.S0_IW*scaleFactor;
fitres.S0_EW                = fitres.S0_EW*scaleFactor;

%%%% Step 6: Save the results
save('MWI_MCR_fitting_result','fitres');
