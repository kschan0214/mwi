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
function fitRes = mwi_3cx_R2s_givenICVFV1_HCFM(algoPara,imgPara)
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
isNormCost = algoPara.isNormCost;

te    = double(imgPara.te);
b0dir = double(imgPara.b0dir);
data  = double(imgPara.img);
mask  = double(imgPara.mask);
fm    = double(imgPara.fieldmap);
pini  = double(imgPara.pini);
icvf  = double(imgPara.icvf);
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
                if mask(ky,kx,kz)>0 && sum(squeeze(ff(ky,kx,kz,:))) > 0
                    % T2*w
                    s       = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0     = fm(ky,kx,kz);
                    pini0   = pini(ky,kx,kz);
                    theta0  = squeeze(theta(ky,kx,kz,:));
                    ff0     = squeeze(ff(ky,kx,kz,:));
                    icvf0   = icvf(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,icvf0,theta0,ff0,db0,pini0,numMagn,isWeighted,isNormCost,userDefine,isInvivo,options,DEBUG);
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
                if mask(ky,kx,kz)>0 && sum(squeeze(ff(ky,kx,kz,:))) > 0
                    % T2*w
                    s       = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                    db0     = fm(ky,kx,kz);
                    pini0   = pini(ky,kx,kz);
                    theta0  = squeeze(theta(ky,kx,kz,:));
                    ff0     = squeeze(ff(ky,kx,kz,:));
                    icvf0   = icvf(ky,kx,kz);
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,icvf0,theta0,ff0,db0,pini0,numMagn,isWeighted,isNormCost,userDefine,isInvivo,options,DEBUG);
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
function [x,res] = FitModel(s,te,icvf0,theta0,ff0,db0,pini0,numMagn,isWeighted,isNormCost,userDefine,isInvivo,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end



% set initial guesses
if isInvivo
    % in vivo reference
    Amy0   = 0.1*abs(s(1));     Amylb   = 0;        Amyub   = 2*abs(s(1));
    Aie0   = 0.9*abs(s(1));  	Aielb   = 0;        Aieub   = 2*abs(s(1));
    r2smw0 = 100;            	r2smwlb = 40;       r2smwub = 300;
    r2siw0 = 16;                r2siwlb = 6;        r2siwub = 40;
    r2sew0 = 21;                r2sewlb = 6;        r2sewub = 40;
    
else
    % ex vivo reference
%     Amy0   = 0.15*abs(s(1));    Amylb   = 0;        Amyub   = 2*abs(s0);
%     Aax0   = 0.65*abs(s(1));  	Aaxlb   = 0;        Aaxub   = 2*abs(s0);
%     Aex0   = 0.2*abs(s(1));    	Aexlb   = 0;        Aexub   = 2*abs(s0);
%     r2smw0 = 10;            	r2smwlb = 1;        r2smwub = 25;
%     r2siw0 = 54;                t2siwlb = 25;       t2siwub = 200;
%     r2sew0 = 38;                t2sewlb = 25;       t2sewub = 200;
end
    
% set initial guess and fitting boundaries
x0 = double([Amy0 ,Aie0 ,r2smw0 ,r2siw0 ,r2sew0 ]);
lb = double([Amylb,Aielb,r2smwlb,r2siwlb,r2sewlb]);
ub = double([Amyub,Aieub,r2smwub,r2siwub,r2sewub]);
if numMagn~=numel(te) % non magnitude fitting

    fbg0   = db0*2*pi;          	fbglb     = (db0-25)*2*pi;     fbgub  = (db0+25)*2*pi;
    if isnan(pini0)
        pini0  = angle(exp(1i*(-2*pi*db0*te(1)-angle(s(1)))));
    end
    pinilb = -2*pi;         piniub=2*pi;

    % extra parameters
    x0 = double([x0,fbg0,pini0]);
    lb = double([lb,fbglb,pinilb]);
    ub = double([ub,fbgub,piniub]);
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
[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,icvf0,theta0,ff0,numMagn,isWeighted,isNormCost,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,icvf,theta,ff,numMagn,isWeighted,isNormCost,DEBUG)
% distribute fitting parameters
Amw=x(1);    Aie=x(2);
t2smw=1/x(3);  t2siw=1/x(4); t2sew=1/x(5);

if numMagn==numel(te) % magnitude fitting
    fbg=0;        pini=0;
else    % other fittings
    fbg=x(6)/(2*pi);     pini=x(7);
end

Aiw = Aie*icvf; Aew = Aie*(1-icvf);

b0    	= 3;
rho_mw  = 0.43;
E       = 0.02;
x_a     = -0.1;
x_i     = -0.1;

sHat = zeros([length(s) length(theta)]);
for kfo = 1:length(theta)
    g   = sqrt(Aiw/(Aiw+Amw/rho_mw));
    [freq_mw,freq_iw] = HCFM_FreqShift(theta(kfo),x_i,x_a,g,E,b0);
    fmwbg = freq_mw + fbg;
    fiwbg = freq_iw + fbg;
%     sin2theta = sin(theta(kfo)).^2;
%     sHat(:,kfo) = mwi_model_ssSPGR_3T2scc_HCFM(te,A_mw,A_iw,A_ew,t2s_mw,t2_fw,sin2theta,param);
    sHat(:,kfo) = mwi_model_3cc_nam2015(te,Amw,Aiw,Aew,t2smw,t2siw,t2sew,fmwbg,fiwbg,fbg,pini);
end
sHat = sum(bsxfun(@times,sHat,ff(:).'),2);


% simulate signal based on parameter input


if size(sHat,1) ~= size(s,1)
    sHat = sHat.';
end

% compute fitting residual
if isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
%     w = sqrt(abs(s)/norm(abs(s(:))));
%     w = sqrt(abs(s)/sum(abs(s(:))));
    w = sqrt(abs(s));
else
    % compute the cost without weights (=same weights)
    w = sqrt(ones(size(s))/numel(s));
end
err = computeFiter(s,sHat,numMagn,w);

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
if isNormCost
    err = err ./ norm(abs(s(:)));
    % err = err ./ mean(abs(s(:)));
end

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
% check # of phase-corrupted echoes
try algoPara2.isNormCost = algoPara.isNormCost;         catch; algoPara2.isNormCost = true; end
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
        error('Fibre orienation map or theta map is required');
    end
end
% b0dir
try
    imgPara2.b0dir = imgPara.b0dir;
    disp('B0dir input: True');
catch
    imgPara2.b0dir = [0,0,1];
    disp('B0dir input: false');
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
if algoPara2.isNormCost
    disp('Cost function is normalised by signal intensity: True');
else
    disp('Cost function is normalised by signal intensity: False');
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