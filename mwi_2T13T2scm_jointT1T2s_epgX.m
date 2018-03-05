%% fitRes = mwi_2T13T2scm_jointT1T2s_epgX(algoPara,imgPara)
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 January 2018
% Date last modified:
%
%
function fitRes = mwi_2T13T2scm_jointT1T2s_epgX(algoPara,imgPara)

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

% capture all fitting settings
maxIter    = algoPara.maxIter;
userDefine  = algoPara.userDefine;
DEBUG      = algoPara.DEBUG;
% isParallel = algoPara.isParallel;

% capture all images related data
te    = imgPara.te;
tr    = imgPara.tr;
fa    = imgPara.fa;
data  = imgPara.img;
mask  = imgPara.mask;
b1map = imgPara.b1map;

[ny,nx,nz,~,~] = size(data);

% Magnitude fitting
numMagn    = length(imgPara.te);

% lsqnonlin setting
fdss = [1e-4,1e-4,1e-4,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8,1e-8];
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'FiniteDifferenceStepSize',fdss,'MaxFunctionEvaluations',200*11,...
        'FunctionTolerance',1e-6,'StepTolerance',1e-1,'Display','off');

% if DEBUG then display fitting message
if DEBUG
    options.Display = 'final';
end

% we are fitting 11 parameters with this model
estimates = zeros(ny,nx,nz,11);
resnorm   = zeros(ny,nx,nz);
for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            if mask(ky,kx,kz)>0
                % 1st dim: T1w; 2nd dim: T2*w
                s = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                b1 = b1map(ky,kx,kz);
                [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,numMagn,userDefine,options,DEBUG);
            end
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b1,numMagn,userDefine,options,DEBUG)
% define initial guesses
% estimate rho0 of the first echo
[~,rho0] = DESPOT1(abs(s(:,1)),fa,tr,'b1',b1);
% estimate t1 from later echo
[t10,~] = DESPOT1(abs(s(:,end-3)),fa,tr,'b1',b1);

% in case rho0 and t10 estimation go wrong then use defualt values
if rho0<0
    rho0=max(abs(s));
end
if t10<0
    t10=1000e-3;
end

Amy0   = 0.1*rho0;            Amylb   = 0;        Amyub   = 1*rho0;
Aax0   = 0.6*rho0;            Aaxlb   = 0;        Aaxub   = 1*rho0;
Aex0   = 0.3*rho0;            Aexlb   = 0;        Aexub   = 1*rho0;
t2smy0 = 10e-3;     t2smylb = 3e-3;     t2smyub = 25e-3;
t2sax0 = 64e-3; 	t2saxlb = 25e-3;    t2saxub = 200e-3;
t2sex0 = 48e-3; 	t2sexlb = 25e-3;    t2sexub = 200e-3;
fmy0   = 5;                fmylb   = 5-75;    fmyub   = 5+75;
fax0   = 0;                 faxlb   = -25;      faxub   = +25;
t1my0  = 300e-3;  	t1mylb  = 50e-3;    t1myub  = 650e-3;
t1l0   = t10;     	t1llb  = 500e-3;    t1lub  = 2000e-3;
kx0 = 2;    kxlb=0;   kxub = 6;

% set initial guess and fitting bounds here
if isempty(userDefine.x0)
    x0 = [Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,t1my0,t1l0,fmy0,fax0,kx0];
else
    x0 = userDefine.x0;
end
if isempty(userDefine.lb)
    lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,t1mylb,t1llb,fmylb,faxlb,kxlb];
else
    lb = userDefine.lb;
end
if isempty(userDefine.ub)
    ub = [Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,t1myub,t1lub,fmyub,faxub,kxub];
else
    ub = userDefine.ub;
end

if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
    x0
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,numMagn,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,numMagn,DEBUG)
% capture all fitting parameters
Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5); t2sex=x(6);
t1my=x(7);  t1l=x(8);
fmy=x(9);  fax=x(10);
kx=x(11);


% sHat = mwi_model_2cm_jointT1T2s_epgX(fa,te,tr,A0,mwf,t2ss,t2sl,t1s,t1l,fs,0,b1,kx);
sHat = mwi_model_2T13T2scm_jointT1T2s_epgX(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1l,fmy,fax,b1,kx);

err = computeFiter(s,sHat,numMagn);

% error weighted by echo strength
% w = abs(s(:))./norm(abs(s(:)));
% err = w.*err;

err = err./norm(abs(s));
% if DEBUG then plots current fitting result
if DEBUG
    global DEBUG_resnormAll
    figure(99);subplot(211);plot(te(:).',abs(permute(s,[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))+10]);
    title('Magnitude');
    plot(te(:).',abs(permute(sHat,[2 1])),'x-');plot(te(:).',(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'ro-.');
    hold off;
    text(te(1)/3,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
    text(te(1)/3,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f,T1my=%f,T1l=%f,kmy=%f',...
        Amy,Aax,Aex,t2smy,t2sax,t2sex,fmy,fax,t1my,t1l,kx));
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

imgPara2 = imgPara;
algoPara2 = algoPara;
isValid = true;

% check if the number of flip angles matches with the data
if length(imgPara.fa) ~= size(imgPara.img,5)
    isValid = false;
end
% check if the number of echo times matches with the data
if length(imgPara.te) ~= size(imgPara.img,4)
    isValid = false;
end

% check maximum iterations allowed
try
    algoPara2.maxIter = algoPara.maxIter;
catch
    algoPara2.maxIter = 500;
end
% check parallel computing 
try
    algoPara2.isParallel = algoPara.isParallel;
catch
    algoPara2.isParallel = false;
end
% check user bounds and initial guesses
try
    algoPara2.userDefine.x0 = algoPara.userDefine.x0;
catch
    algoPara2.userDefine.x0 = [];
end
try
    algoPara2.userDefine.lb = algoPara.userDefine.lb;
catch
    algoPara2.userDefine.lb = [];
end
try
    algoPara2.userDefine.ub = algoPara.userDefine.ub;
catch
    algoPara2.userDefine.ub = [];
end

% check signal mask
try
    imgPara2.mask = imgPara.mask;
catch
    imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
end

% check b1 map
try
    imgPara2.b1map = imgPara.b1map;
catch
    imgPara2.b1map = ones(size(imgPara.img,1),size(imgPara.img,2),size(imgPara.img,3));
end

end