%% fitRes = mwi_3cm_jointT1T2s_2T1(algoPara,imgPara)
%
% Input
% --------------
% algoPara.maxIter    : maximum no. of iterations allowed (default: 500)
% algoPara.isParallel : parallel computing(parfor) (default: false)
% algoPara.usrDefine  : user defined fitting initial guesses and bounds
% imgPara.te          : echo times (s)
% imgPara.tr          : repetition times (s)
% imgPara.fa          : flip angles (degree)
% imgPara.img         : 5D images, [ny,ny,nz,necho,nfa]
% imgPara.mask        : signal mask (default: rel. signal>0.05)
% imgPara.b1map       : B1 map (default: 1)
%
% Output
% --------------
% fitRes.estimates    : estimates,[Amy,Aax,Aex,T2smy,T2sax,T2sex,T1my,T1l,fmy,fax]
% fitRes.resnorm      : fitting residual norm
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 January 2018
% Date last modified: 21 February 2018
%
%
function fitRes = mwi_3cm_jointT1T2s_2T1(algoPara,imgPara)

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

% capture all fitting settings
maxIter    = algoPara.maxIter;
usrDefine  = algoPara.usrDefine;
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
fdss = [1e-4,1e-4,1e-4,1e-7,1e-7,1e-7,1e-6,1e-6,1e-5,1e-5];
options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'FiniteDifferenceStepSize',fdss,'MaxFunctionEvaluations',200*11);

% if DEBUG then display fitting message
if DEBUG
    options.Display = 'off';
end

% we are fitting 10 parameters with this model
estimates = zeros(ny,nx,nz,10);
resnorm   = zeros(ny,nx,nz);
for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            if mask(ky,kx,kz)>0
                % 1st dim: T1w; 2nd dim: T2*w
                s = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                b1 = b1map(ky,kx,kz);
                [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,numMagn,usrDefine,options,DEBUG);
            end
        end
    end
end

% store fitting results into structure format
fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b1,numMagn,usrDefine,options,DEBUG)
% if DEBUG then create an array to store resnorm of all iterations
if DEBUG
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end
% define initial guesses
% estimate m0 of the first echo
[~,m0] = DESPOT1(abs(s(:,1)),fa,tr,'b1',b1);
% estimate t1 from later echo
[t10,~] = DESPOT1(abs(s(:,end-3)),fa,tr,'b1',b1);

% in case m0 and t10 estimation go wrong then use defualt values
if m0<0
    m0=max(abs(s));
end
if t10<0
    t10=1000e-3;
end

Amy0   = 0.1*m0;	Amylb   = 0;        Amyub   = 0.3*m0;
Aax0   = 0.6*m0;	Aaxlb   = 0;        Aaxub   = 1*m0;
Aex0   = 0.3*m0;	Aexlb   = 0;        Aexub   = 1*m0;
t2smy0 = 10e-3;     t2smylb = 3e-3;     t2smyub = 25e-3;
t2sax0 = 64e-3; 	t2saxlb = 25e-3;    t2saxub = 500e-3;
t2sex0 = 48e-3; 	t2sexlb = 25e-3;    t2sexub = 500e-3;
t1my0  = 200e-3;  	t1mylb  = 50e-3;    t1myub  = 650e-3;
t1l0   = t10;     	t1llb  = 650e-3;    t1lub  = 2000e-3;
fmy0   = 5;      	fmylb   = 5-75;     fmyub   = 5+75;
fax0   = 0;       	faxlb   = -25;      faxub   = 25;

% set initial guess and fitting bounds here
if isempty(usrDefine.x0)
    x0 = [Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,t1my0,t1l0,fmy0,fax0];
else
    x0 = usrDefine.x0;
end
if isempty(usrDefine.lb)
    lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,t1mylb,t1llb,fmylb,faxlb];
else
    lb = usrDefine.lb;
end
if isempty(usrDefine.ub)
    ub = [Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,t1myub,t1lub,fmyub,faxub];
else
    ub = usrDefine.ub;
end

% run lsqnonlin!
[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,numMagn,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,numMagn,DEBUG)
% capture all fitting parameters
Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5);  t2sex=x(6);
t1my=x(7);  t1l=x(8);
fmy=x(9);  fax=x(10);

% simulate signal based on input parameters
sHat = mwi_model_3cm_jointT1T2s_2T1(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1l,fmy,fax,b1);

% compute the residual between the simulated signal and measured signal
err = computeFiter(s,sHat,numMagn);

% if DEBUG then plots current fitting result
if DEBUG
    global DEBUG_resnormAll
%     w = abs(s(:))/norm(abs(s(:)));
%         w = ones(size(s));
%         w(1,4,:) = 100;
%         w = bsxfun(@times,abs(s),1./vecnorm(abs(s),2,2));
%         w = w(:)./norm(w(:));
%     err = err.*w(:);
    figure(99);subplot(211);plot(te(:).',abs(permute(s,[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))+10]);
    title('Magnitude');
    plot(te(:).',abs(permute(sHat,[2 1])),'x-');plot(te(:).',(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'ro-.');
    hold off;
    text(te(1)*1.1,max(abs(s(:))*0.2),sprintf('resnorm=%f',sum(err(:).^2)));
    text(te(1)*1.1,max(abs(s(:))*0.1),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f,T1my=%f,T1l=%f',...
        Amy,Aax,Aex,t2smy,t2sax,t2sex,fmy,fax,t1my,t1l));
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

% check DEBUG
try
    algoPara2.DEBUG = algoPara.DEBUG;
catch
    algoPara2.DEBUG = false;
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
    algoPara2.usrDefine.x0 = algoPara.usrDefine.x0;
catch
    algoPara2.usrDefine.x0 = [];
end
try
    algoPara2.usrDefine.lb = algoPara.usrDefine.lb;
catch
    algoPara2.usrDefine.lb = [];
end
try
    algoPara2.usrDefine.ub = algoPara.usrDefine.ub;
catch
    algoPara2.usrDefine.ub = [];
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