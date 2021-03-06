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
% Date last modified: 28 January 2019
%
%
function fitRes = mwi_3T13T2s_cx_2steps_linearX(algoPara,imgPara)
% make sure computeFiter in the same tool is used
% addpath('utils/');

disp('Myelin water imaing: VFA-T2* model');
% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

% get debug mode and verbose
DEBUG   = algoPara.DEBUG;
verbose = algoPara.verbose;

% capture all fitting settings
maxIter      = algoPara.maxIter;
fcnTol       = algoPara.fcnTol;
numMagn      = algoPara.numMagn;
isWeighted   = algoPara.isWeighted;
weightMethod = algoPara.weightMethod;
userDefine   = algoPara.userDefine;
isParallel   = algoPara.isParallel;
isInvivo     = algoPara.isInvivo;

% check the user weighting method matches the available one or not
if strcmpi(weightMethod,'norm') && strcmpi(weightMethod,'1stEcho')
    disp('Do not support the input weighting method');
    weightMethod = 'norm';
end

% if DEBUG on then disables parallel computing
if DEBUG
    isParallel = false;
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

%%%%%%%%%% display fitting message %%%%%%%%%%
if verbose
    disp('The following fitting parameters are used:');
    fprintf('Max. iterations = %i\n',maxIter);
    fprintf('Function tolerance = %e\n',fcnTol);
   
    
    if isWeighted(1)
        disp('MC-T2* fitting: Weighted cost function: True');
        disp(['Weighting method: ' weightMethod]);
    else
        disp('MC-T2* fitting: Weighted cost function: False');
    end
    if isWeighted(2)
        disp('EPG-X fitting: Weighted cost function: True');
        disp(['Weighting method: ' weightMethod]);
    else
        disp('EPG-X fitting: Weighted cost function: False');
    end
    
    if isInvivo
        disp('Initial guesses for in vivo study');
    else
        disp('Initial guesses for ex vivo study');
    end
    
    % type of fitting
    if numMagn(1)==0
        disp('MC-T2* fitting: Fitting complex model with complex data');
    elseif numMagn(1)==numel(te)
        disp('MC-T2* fitting: Fitting complex model with magnitude data');
    else
        fprintf('MC-T2* fitting: Fitting complex model with %i magnitude data and %i complex data\n',numMagn,numel(te)-numMagn);
    end
    % type of fitting
    if numMagn(2)==0
        disp('EPG-X fitting: Fitting complex model with complex data');
    elseif numMagn(2)==numel(te)
        disp('EPG-X fitting: Fitting complex model with magnitude data');
    else
        fprintf('EPG-X fitting: Fitting complex model with %i magnitude data and %i complex data\n',numMagn,numel(te)-numMagn);
    end
    
    % initial guess and fitting bounds
    if isempty(userDefine(1).x0)
        disp('MC-T2* fitting: Default initial guess: True');
    else
        disp('MC-T2* fitting: Default initial guess: False');
    end
    if isempty(userDefine(1).lb)
        disp('MC-T2* fitting: Default lower bound: True');
    else
        disp('MC-T2* fitting: Default lower bound: False');
    end
    if isempty(userDefine(1).ub)
        disp('MC-T2* fitting: Default upper bound: True');
    else
        disp('MC-T2* fitting: Default upper bound: False');
    end
    if isempty(userDefine(2).x0)
        disp('EPG-X fitting: Default initial guess: True');
    else
        disp('EPG-X fitting: Default initial guess: False');
    end
    if isempty(userDefine(2).lb)
        disp('EPG-X fitting: Default lower bound: True');
    else
        disp('EPG-X fitting: Default lower bound: False');
    end
    if isempty(userDefine(2).ub)
        disp('EPG-X fitting: Default upper bound: True');
    else
        disp('EPG-X fitting: Default upper bound: False');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny,nx,nz,~,nfa] = size(data);

% lsqnonlin setting
options_step1 = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*11,...
        'FunctionTolerance',fcnTol(1),'StepTolerance',1e-6);
options_step2 = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*11,...
        'FunctionTolerance',fcnTol(2),'StepTolerance',1e-6);
  
% fdss = [1e-4,1e-4,1e-4,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8];  
% options.FiniteDifferenceStepSize = fdss;

% if DEBUG then display fitting message
if ~DEBUG
    options_step1.Display = 'off';
    options_step2.Display = 'off';
end


% EPG-X has an extra exchange term
if numMagn(2)==numel(te)   % magnitude fitting has 11 estimates
    estimates_step2 = zeros(ny,nx,nz,8);
else
    estimates_step2 = zeros(ny,nx,nz,8+nfa*2); % others have 12 estimates
end


if numMagn(1)==numel(te)   % magnitude fitting has 10 estimates
    estimates_step1 = zeros(ny,nx,nz,nfa*3+5);
else
    estimates_step1 = zeros(ny,nx,nz,nfa*5+5);
end

resnorm_step1 = zeros(ny,nx,nz);
resnorm_step2 = zeros(ny,nx,nz);
if isParallel
    for kz=1:nz
        if verbose
            fprintf('Processing slice %i...\n',kz);
        end
%         numVoxelToBeFitted = length(find(mask(:,:,kz)==1));
%         fprintf('%i voxles need to be fitted...\n',numVoxelToBeFitted);
%         fprintf('Progress (%%): ');
        for ky=1:ny
            parfor kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s       = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1      = b1map(ky,kx,kz);
                    db0     = squeeze(fm(ky,kx,kz,:));
                    pini0   = squeeze(pini(ky,kx,kz,:));
                    [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,db0,pini0,npulse,numMagn,isWeighted,weightMethod,isInvivo,model,userDefine,options,DEBUG);
                end
            end
%             if numVoxelToBeFitted > 0
%                 numVoxelFitted = length(find(mask(1:ky,:,kz)==1));
%                 percentFinish = floor(numVoxelFitted/numVoxelToBeFitted * 100);
%                 if percentFinish>1
%                   for j=0:log10(percentFinish-1)
%                       fprintf('\b'); % delete previous counter display
%                   end
%                   fprintf('%% %i ',percentFinish);
%                 end
%             end
            if mod(ky,5) == 0;
                fprintf('%i ', ky);
            end
        end
        fprintf('\n');
    end
else
    for kz=1:nz
        if verbose
            fprintf('Processing slice %i\n',kz);
        end
        for ky=1:ny
            for kx=1:nx
                if mask(ky,kx,kz)>0
                    % 1st dim: T1w; 2nd dim: T2*w
                    s = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                    b1 = b1map(ky,kx,kz);
                    db0     = squeeze(fm(ky,kx,kz,:));
                    pini0   = squeeze(pini(ky,kx,kz,:));
                    [estimates_step1(ky,kx,kz,:),resnorm_step1(ky,kx,kz)] = FitModel_t2s(s,te,db0,pini0,numMagn(1),isWeighted(1),weightMethod,isInvivo,userDefine(1),options_step1,DEBUG);
                    [estimates_step2(ky,kx,kz,:),resnorm_step2(ky,kx,kz)] = FitModel_vfa(s,fa,te,tr,b1,squeeze(estimates_step1(ky,kx,kz,:)),numMagn(2),isWeighted(2),weightMethod,isInvivo,userDefine(2),options_step2,DEBUG);
                end
            end
            if mod(ky,5) == 0;
                fprintf('%i ', ky);
            end
        end
        fprintf('\n');
    end
end

fitRes.estimates_step1 = estimates_step1;
fitRes.resnorm_step1   = resnorm_step1;
fitRes.estimates_step2 = estimates_step2;
fitRes.resnorm_step2   = resnorm_step2;
fitRes.estimates       = cat(4,estimates_step2(:,:,:,1:4),...
                               estimates_step1(:,:,:,3*nfa+1:3*nfa+3),...
                               estimates_step2(:,:,:,5:6),...
                               estimates_step1(:,:,:,3*nfa+4:3*nfa+5),...
                               estimates_step2(:,:,:,7:end));

end

%% Setup lsqnonlin and fit with the MC-T2s model
function [x,res] = FitModel_t2s(s,te,db0,pini,numMagn,isWeighted,weightMethod,isInvivo,userDefine,options,DEBUG)
% define initial guesses
% quick estimate T1w signal for each flip angle
s0 = zeros(1,size(s,1));
for kfa=1:size(s,1)
    b = [ones(length(te),1), -te(:)]\log(abs(s(kfa,:).'));
    % r2s(kx,ky,kz) = b(2);
    s0(kfa) = exp(b(1));
    if s0(kfa)<0
        s0(kfa)=0;
    end
end

if isInvivo
    % common initial guesses for in vivo study
    Amw0   = 0.1*abs(s(:,1)).';	Amwlb   = zeros(1,size(s,1));	Amwub   = 2*s0;
    Aiw0   = 0.6*abs(s(:,1)).';	Aiwlb   = zeros(1,size(s,1));	Aiwub   = 2*s0;
    Aew0   = 0.3*abs(s(:,1)).';	Aewlb   = zeros(1,size(s,1));	Aewub   = 2*s0;
    t2smw0 = 10e-3;            	t2smwlb = 1e-3;               	t2smwub = 25e-3;
    t2siw0 = 64e-3;            	t2siwlb = 25e-3;                t2siwub = 200e-3;
    t2sew0 = 192e-3;           	t2sewlb = 0e-3;                t2sewub = 400e-3;
    
else
    % common initial guesses for ex vivo study
    Amw0   = 0.15*abs(s(:,1)).';	Amwlb   = zeros(1,size(s,1));	Amwub   = 2*s0;
    Aiw0   = 0.6*abs(s(:,1)).';     Aiwlb   = zeros(1,size(s,1)); 	Aiwub   = 2*s0;
    Aew0   = 0.25*abs(s(:,1)).';    Aewlb   = zeros(1,size(s,1)); 	Aewub   = 2*s0;
    t2smw0 = 10e-3;                 t2smwlb = 1e-3;                	t2smwub = 20e-3;
    t2siw0 = 54e-3;                 t2siwlb = 20e-3;                t2siwub = 200e-3;
    t2sew0 = 38e-3;                 t2sewlb = 20e-3;                t2sewub = 200e-3;
    
end

fmw0   = 5;  	fmwlb   = -75;      fmwub   = 75;
fiw0   = 0;  	fiwlb   = -25;  	fiwub   = 25;

if numMagn==numel(te) % magnitude fitting
    
    % set initial guess and fitting boundaries
    x0 = double([Amw0 ,Aiw0 ,Aew0 ,t2smw0 ,t2siw0 ,t2sew0 ,fmw0 ,fiw0 ]);
    lb = double([Amwlb,Aiwlb,Aewlb,t2smwlb,t2siwlb,t2sewlb,fmwlb,fiwlb]);
    ub = double([Amwub,Aiwub,Aewub,t2smwub,t2siwub,t2sewub,fmwub,fiwub]);
    
else    % other fittings
    
    % total field in Hz, initial phase in radian
    totalField0 = db0(:).';  totalFieldlb   = db0(:).'-100;                 totalFieldub    = db0(:).'+100;
    pini0       = pini(:).'; pinilb         = ones(size(pini0))*(-2*pi);    piniub          = ones(size(pini0))*2*pi;
    
    % set initial guess and fitting boundaries
    x0 = double([Amw0 ,Aiw0 ,Aew0 ,t2smw0 ,t2siw0 ,t2sew0 ,fmw0 ,fiw0 ,totalField0,pini0]);
    lb = double([Amwlb,Aiwlb,Aewlb,t2smwlb,t2siwlb,t2sewlb,fmwlb,fiwlb,totalFieldlb,pinilb]);
    ub = double([Amwub,Aiwub,Aewub,t2smwub,t2siwub,t2sewub,fmwub,fiwub,totalFieldub,piniub]);
    
end

% set user defined initial guess and fitting bounds here
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
    global DEBUG_resnormAll x1
    x1 = [];
    DEBUG_resnormAll=[];
    x0
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc_t2s(y,s,te,numMagn,isWeighted,weightMethod,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc_t2s(x,s,te,numMagn,isWeighted,weightMethod,DEBUG)
nfa = size(s,1);
% capture all fitting parameters
Amw=x(1:nfa);       Aiw=x(nfa+1:2*nfa);	Aew=x(2*nfa+1:3*nfa);
t2smw=x(3*nfa+1);   t2iw=x(3*nfa+2);   t2_prime=x(3*nfa+3);
fmw=x(3*nfa+4);     fiw=x(3*nfa+5); 

t2sew = 1/ (1/t2iw+1/t2_prime);

if numMagn==numel(te) % magnitude fitting      
    % no initial phase 
    pini        = zeros(1,nfa);
    totalfield  = zeros(1,nfa);
else    % other fittings    
    totalfield = x(3*nfa+6:3*nfa+6+nfa-1);
    pini       = x(3*nfa+6+nfa:end);
end

% simulate signal based on parameter input
sHat = zeros(size(s));
for kfa = 1:nfa
    sHat(kfa,:) = mwi_model_3cc_nam2015(te,Amw(kfa),Aiw(kfa),Aew(kfa),...
                                           t2smw,t2iw,t2sew,...
                                           fmw+totalfield(kfa),fiw+totalfield(kfa),totalfield(kfa),pini(kfa));
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
    w = sqrt(ones(size(s))/numel(s));
    % compute the cost without weights
    err = computeFiter(s,sHat,numMagn,w);
end

% residual normalied by measured signal
% 20180316 TODO:maybe for weighted cost there should be another way to do this 
err = err./norm(abs(s));

% if DEBUG then plots current fitting result
if 0
    global DEBUG_resnormAll x1
    x1 = [x1,x(:)];
    figure(99);subplot(211);plot(te(:).',abs(permute(s,[2 1])),'^-');hold on;ylim([0-min(abs(s(:))),max(abs(s(:)))+10]);
    title('Magnitude');
    plot(te(:).',abs(permute(sHat,[2 1])),'x-.');plot(te(:).',(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'o-.');
    hold off;
    text(te(1)/3,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
    text(te(1)/3,max(abs(s(:))*0.1),sprintf('t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f',...
        t2smw,t2iw,t2sew,fmw,fiw));
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    subplot(212);plot(DEBUG_resnormAll);
    if length(DEBUG_resnormAll) <300
        xlim([0 300]);
    else
        xlim([length(DEBUG_resnormAll)-300 length(DEBUG_resnormAll)]);
    end
    drawnow;
end

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel_vfa(s,fa,te,tr,b1,estimates,numMagn,isWeighted,weightMethod,isInvivo,userDefine,options,DEBUG)
nfa = length(fa);
% get T1w signal estimated from MC-T2s fitting
Amw = estimates(1:nfa);
Aiw = estimates(nfa+1:2*nfa);
Aew = estimates(2*nfa+1:3*nfa);

% T2* and frequencies are fixed in this step
t2s = estimates(3*nfa+1:3*nfa+3);
freq = estimates(3*nfa+4:3*nfa+5);
if length(estimates) < 3*nfa+6
    totalField = zeros(1,nfa);
    pini = zeros(1,nfa);
else
    totalField = estimates(3*nfa+6:4*nfa+5).';	
    pini       = estimates(4*nfa+6:end).';
end

% estimate proton density of each pool based on MC-T2s fitting
rho = zeros(1,3);
t1 = zeros(1,3);
[t1(1),rho(1)] = DESPOT1(abs(Amw),fa,tr,'b1',b1);
[t1(2),rho(2)] = DESPOT1(abs(Aiw),fa,tr,'b1',b1);
[t1(3),rho(3)] = DESPOT1(abs(Aew),fa,tr,'b1',b1);

% estimate proton density of total water using 1st echo
[~,rho0] = DESPOT1(abs(s(:,1)),fa,tr,'b1',b1);
% estimate long t1 from later echoes
[t10,~] = DESPOT1(abs(s(:,end-3)),fa,tr,'b1',b1);

% in case rho0 and t10 go wrong then use defualt values
if rho0<0
    rho0=max(abs(s(:)));
end
if t10<0 || t10>3000e-3
    t10=1000e-3;
end

Amw0   = rho(1); 	Amwlb   = 0;        Amwub   = 2*rho0;
Aiw0   = rho(2); 	Aiwlb   = 0;        Aiwub   = 2*rho0;
Aew0   = rho(3); 	Aewlb   = 0;        Aewub   = 2*rho0;
% Amw0   = 0.2*rho0; 	Amwlb   = 0;        Amwub   = 2*rho0;
% Aiw0   = 0.65*rho0; 	Aiwlb   = 0;        Aiwub   = 2*rho0;
% Aew0   = 0.35*rho0; 	Aewlb   = 0;        Aewub   = 2*rho0;
% Ass0   = 0.24*sum(rho0/(1-0.24)); 	Asslb   = 0;        Assub   = 2*rho0;
if isInvivo
    % common initial guesses for in vivo study
    t1mw0   = 228e-3;  	t1mwlb   = 25e-3;  	t1mwub   = 650e-3;
    t1iw0   = 1.5;     	t1iwlb   = 500e-3; 	t1iwub   = 4000e-3;
    t1ew0   = 1.5;      t1ewlb  = 500e-3;   t1ewub  = 4;
    kba0    = 0.6;       	kbalb    = 0;     kbaub    = 20;       % exchange rate from long T1 to short T1
    kbc0    = 6.7;       	kbclb    = 0;      	kbcub    = 20;       % exchange rate from long T1 to short T1

else
    % common initial guesses for ex vivo study
    t1s0   = 100e-3;  	t1slb   = 50e-3;  	t1sub   = 400e-3;
    t1iw0   = t10;     	t1iwlb   = 300e-3; 	t1iwub   = 2000e-3;
    kba0   = 2;       	kowmwlb   = 0;      	kowmwub   = 6;       % exchange rate from long T1 to short T1
    
end

% set initial guess and fitting boundaries
% x0 = double([Amw0 ,Aiw0 ,Aew0 ,Ass0,t1s0 ,t1l0,kowmw0,kmwss0]);
% lb = double([Amwlb,Aiwlb,Aewlb,Asslb,t1slb,t1llb,kowmwlb,kmwsslb]);
% ub = double([Amwub,Aiwub,Aewub,Assub,t1sub,t1lub,kowmwub,kmwssub]);
x0 = double([Amw0 ,Aiw0 ,Aew0 ,t1mw0,t1iw0,t1ew0,kba0,kbc0]);
lb = double([Amwlb,Aiwlb,Aewlb,t1mwlb,t1iwlb,t1ewlb,kbalb,kbclb]);
ub = double([Amwub,Aiwub,Aewub,t1mwub,t1iwub,t1ewub,kbaub,kbcub]);
% x0 = double([Amw0 ,Aiw0 ,Aew0  ,t1iw0,kba0]);
% lb = double([Amwlb,Aiwlb,Aewlb,t1iwlb,kowmwlb]);
% ub = double([Amwub,Aiwub,Aewub,t1iwub,kowmwub]);

if numMagn==numel(te) % magnitude fitting

    
else    % other fittings
    totalField0 = totalField(:).';	totalFieldlb   = totalField(:).'-100;	totalFieldub    = totalField(:).'+100;
    pini0       = pini(:).';   	pinilb         = ones(size(pini0))*(-2*pi);         piniub          = ones(size(pini0))*2*pi;

        % set initial guess and fitting boundaries
        x0 = double([x0 ,totalField0,pini0]);
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
    global DEBUG_resnormAll x2
    DEBUG_resnormAll=[];x2=[];
    x0
end


% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc_vfa(y,s,fa,te,tr,b1,t2s,freq,numMagn,isWeighted,weightMethod,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc_vfa(x,s,fa,te,tr,b1,t2s,freq,numMagn,isWeighted,weightMethod,DEBUG)
% fixed parameters
t2smw=t2s(1); t2siw=t2s(2); t2sew=1/ (1/t2siw+1/t2s(3));
freqmw=freq(1);   freqiw=freq(2); 

% fitting parameters
Amw=x(1);   Aiw=x(2);   Aew=x(3);
t1mw=x(4);
t1iw=x(5);
t1ew=x(6);


kba=x(7);  
kbc=x(8);

% t1ss=228e-3; 
if numMagn==numel(te) % magnitude fitting      
    % no initial phase 
    pini=zeros(size(fa));
    fbg = zeros(size(fa));
    
else    % other fittings
        fbg = x(9:9+length(fa)-1);
        pini       = x(9+length(fa):end);
end

% simulate signal based on parameter input
fa = fa*b1;

% linear system
kac = 0;
kca = 0;
[~,ssHat] = mwi_model_ssSPGR_3T1(fa,tr,...
                              Aiw,Amw,Aew,...
                              t1iw,t1mw,t1ew,...
                              kba,kbc,kca,kac);


% [~,ssHat] = mwi_model_ssSPGR_2T1(fa,tr,...
%                                   Amw,Aew+Aiw,...
%                                   t1ss,t1w,...
%                                   kba);
                          
sHat = zeros(size(s));
for kfa = 1:length(fa)
    Amw_curr = ssHat(2,kfa); fmwbg = freqmw + fbg(kfa);
    Aiw_curr = ssHat(1,kfa); fiwbg = freqiw + fbg(kfa);
    Aew_curr = ssHat(3,kfa); 
   
    
    sHat(kfa,:) = mwi_model_3cc_nam2015(te,Amw_curr,Aiw_curr,Aew_curr,t2smw,t2siw,t2sew,fmwbg,fiwbg,fbg(kfa),pini(kfa));
end

% sHat = mwi_model_3T13T2scc_linearX(fa,te,tr,...
%                                       Amw,Aiw,Aew,Ass,...
%                                       t2smw,t2siw,t2sew,...
%                                       t1w,t1w,t1ss,...
%                                       kowmw,kmwss,...
%                                       freqmw,freqiw,...
%                                       totalfield,pini,b1);

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
    w = sqrt(ones(size(s))/numel(s));
    % compute the cost without weights
    err = computeFiter(s,sHat,numMagn,w);
end

% residual normalied by measured signal
% 20180316 TODO:maybe for weighted cost there should be another way to do this 
err = err./norm(abs(s));

% if DEBUG then plots current fitting result
if DEBUG
    global DEBUG_resnormAll x2
    x2 = [x2,x(:)];
    figure(97);subplot(211);plot(te(:).',abs(permute(s,[2 1])),'^-');hold on;ylim([0-min(abs(s(:))),max(abs(s(:)))+10]);
    title('Magnitude');
    plot(te(:).',abs(permute(sHat,[2 1])),'x-.');plot(te(:).',(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'o-.');
    hold off;
    text(te(1)/3,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
    text(te(1)/3,max(abs(s(:))*0.1),sprintf('Amw=%f,Aiw=%f,Aew=%f,T1mss=%f,T1w=%f,kowmw=%f,kmwss=%f',...
        Amw,Aiw,Aew,t1ew,t1iw,kba,t1mw));
    for kfa = 1:length(fa)
        text(te(1)/3,abs(s(kfa,1)),['FA ' num2str(fa(kfa))]);
    end
    DEBUG_resnormAll = [DEBUG_resnormAll;sum(err(:).^2)];
    subplot(212);plot(DEBUG_resnormAll);
    if length(DEBUG_resnormAll) <300
        xlim([0 300]);
    else
        xlim([length(DEBUG_resnormAll)-300 length(DEBUG_resnormAll)]);
    end
    drawnow;
end

end

%% check and set default
function [algoPara2,imgPara2,isValid]=CheckAndSetPara(algoPara,imgPara)
% copy input to output
imgPara2 = imgPara;
algoPara2 = algoPara;
isValid = true;

% check if the number of flip angles matches with the data's 5th dim
if length(imgPara.fa) ~= size(imgPara.img,5)
    isValid = false;
end
% check if the number of echo times matches with the data's 4th dim
if length(imgPara.te) ~= size(imgPara.img,4)
    isValid = false;
end

% check debug
try
    algoPara2.DEBUG = algoPara.DEBUG;
catch
    algoPara2.DEBUG = false;
end
% check verbose
try
    algoPara2.verbose = algoPara.verbose;
catch
    algoPara2.verbose = true;
end

% check maximum iterations allowed
try
    algoPara2.maxIter = algoPara.maxIter;
catch
    algoPara2.maxIter = 500;
end
% check function tolerance
try
    algoPara2.fcnTol = algoPara.fcnTol;
    if numel(algoPara2.fcnTol) == 1
        algoPara2.fcnTol(2) = algoPara2.fcnTol;
    end
catch
    algoPara2.fcnTol(1:2) = 1e-5;
end
% check parallel computing 
try
    algoPara2.isParallel = algoPara.isParallel;
catch
    algoPara2.isParallel = false;
end
% check weighted sum of cost function
try
    algoPara2.weightMethod = algoPara.weightMethod;
catch
    algoPara2.weightMethod = 'norm';
end
% check method for weighting
try
    algoPara2.isWeighted = algoPara.isWeighted;
    if numel(algoPara2.isWeighted) == 1
        algoPara2.isWeighted(2) = algoPara2.isWeighted;
    end
catch
    algoPara2.isWeighted(1:2) = false;
end
% check # of phase-corrupted echoes
try
    algoPara2.numMagn = algoPara.numMagn;
    if numel(algoPara2.numMagn) == 1
        algoPara2.numMagn(2) = algoPara2.numMagn;
    end
catch
    algoPara2.numMagn(1:2) = numel(imgPara.te);
end
% check user bounds and initial guesses
try
    algoPara2.userDefine(1).x0 = algoPara.userDefine(1).x0;
catch
    algoPara2.userDefine(1).x0 = [];
end
try
    algoPara2.userDefine(1).lb = algoPara.userDefine(1).lb;
catch
    algoPara2.userDefine(1).lb = [];
end
try
    algoPara2.userDefine(1).ub = algoPara.userDefine(1).ub;
catch
    algoPara2.userDefine(1).ub = [];
end
try
    algoPara2.userDefine(2).x0 = algoPara.userDefine(2).x0;
catch
    algoPara2.userDefine(2).x0 = [];
end
try
    algoPara2.userDefine(2).lb = algoPara.userDefine(2).lb;
catch
    algoPara2.userDefine(2).lb = [];
end
try
    algoPara2.userDefine(2).ub = algoPara.userDefine(2).ub;
catch
    algoPara2.userDefine(2).ub = [];
end

% check initial guesses
try
    imgPara2.isInvivo = imgPara.isInvivo;
catch
    imgPara2.isInvivo = true;
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

try
    algoPara2.model = algoPara.model;
    if ~strcmpi(algoPara2.model,'epgx') && ~strcmpi(algoPara2.model,'epg') && ~strcmpi(algoPara2.model,'standard')
        error('Your input model is not supported. Please choose either ''epgx'', ''epg'' or ''standard''. ');
    end
catch
    algoPara2.model = 'epgx';
end

try
    algoPara2.RFphi = algoPara.RFphi;
catch
    algoPara2.RFphi = 50;
end


end