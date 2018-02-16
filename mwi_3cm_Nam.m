%% fitRes = mwi_3cm_jointT1T2s(algoPara,imgPara)
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
function fitRes = mwi_3cm_Nam(algoPara,imgPara)

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

DEBUG = algoPara.DEBUG;

% capture all parameters
numMagn    = 0;
maxIter    = algoPara.maxIter;
% isParallel = algoPara.isParallel;

te    = imgPara.te;
data  = imgPara.img;
mask  = imgPara.mask;
fm    = imgPara.fieldmap;

[ny,nx,nz,~,~] = size(data);

if DEBUG
%     fdss = [1e-2,1e-2,1e-3,1e-7,1e-7,1e-7,1e-6,1e-6,1e-6,1e-5,1e-5];
%     fdss = [1e-4,1e-3,1e-3,1e-7,1e-7,1e-7,1e-4,1e-4,1e-4,1e-4];
    fdss = [1e-4,1e-7,1e-4,1e-3,1e-7,1e-4,1e-3,1e-7,1e-4,1e-4];
    options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'FiniteDifferenceStepSize',fdss,'MaxFunctionEvaluations',200*15,...
        'StepTolerance',1e-6,'FunctionTolerance',1e-6);
else
    options = optimoptions(@lsqnonlin,'Display','off','MaxIter',maxIter);
end

estimates = zeros(ny,nx,nz,10);
resnorm   = zeros(ny,nx,nz);
for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            if mask(ky,kx,kz)>0
                % T2*w
                s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                db0 = fm(ky,kx,kz);
                [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,numMagn,options,DEBUG);
            end
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,te,db0,numMagn,options,DEBUG)
[~,~,s0]=R2starmapping_trapezoidal_voxel(abs(s),te);
Amy0   = 0.1*abs(s0);            Amylb   = 0;        Amyub   = 0.3*abs(s0);
Aax0   = 0.6*abs(s0);            Aaxlb   = 0;        Aaxub   = 1*abs(s0);
Aex0   = 0.3*abs(s0);            Aexlb   = 0;        Aexub   = 1*abs(s0);
% Amy0   = 0.15*abs(s(1));            Amylb   = 0;        Amyub   = 0.3*abs(s(1));
% Aax0   = 0.6*abs(s(1));            Aaxlb   = 0;        Aaxub   = 1*abs(s(1));
% Aex0   = 0.3*abs(s(1));            Aexlb   = 0;        Aexub   = 1*abs(s(1));
t2smy0 = 10e-3;                    t2smylb = 3e-3;     t2smyub = 24e-3;
t2sax0 = 64e-3;               t2saxlb = 24e-3;    t2saxub = 350e-3;
t2sex0 = 48e-3;               t2sexlb = 24e-3;    t2sexub = 350e-3;
fmy0   = db0;                 fmylb   = db0-75;    fmyub   = db0+75;
fax0   = db0;                 faxlb   = db0-25;      faxub   = db0+25;
fex0   = db0;                 fexlb   = db0-25;      fexub   = db0+25;
pini0   = 2*pi*db0*te(1)-angle(s(1));        pinilb = -pi;         piniub=pi;


x0 = double([Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,fmy0,fax0,fex0,pini0]);
lb = double([Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,fmylb,faxlb,fexlb,pinilb]);
ub = double([Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,fmyub,faxub,fexub,piniub]);
% x0 = double([Amy0,t2smy0,fmy0,Aex0,t2sex0,fex0,Aax0,t2sax0,fax0,pini0]);
% lb = double([Amylb,t2smylb,fmylb,Aexlb,t2sexlb,fexlb,Aaxlb,t2saxlb,faxlb,pinilb]);
% ub = double([Amyub,t2smyub,fmyub,Aexub,t2sexub,fexub,Aaxub,t2saxub,faxub,piniub]);


[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,numMagn,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,numMagn,DEBUG)
Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5); t2sex=x(6);
fmybg=x(7);  faxbg=x(8);    fexbg=x(9);
pini=x(10);
% Amy=x(1);   Aax=x(7);   Aex=x(4);
% t2smy=x(2); t2sax=x(8); t2sex=x(5);
% fmybg=x(3);  faxbg=x(9);    fexbg=x(6);
% pini=x(10);

sHat = mwi_model_3cc_nam2015(te,Amy,Aax,Aex,t2smy,t2sax,t2sex,fmybg,faxbg,fexbg,pini);

err = computeFiter(s,sHat,numMagn);

if DEBUG
%     % weighted by signal magnitude
%     w = abs(s(:))/norm(abs(s(:)));
%     err = err.*repmat(w(:),2,1);
    figure(99);subplot(411);plot(te(:).',real(permute(s,[2 1])),'k^-');hold on;ylim([min(real(s(:)))-10,max(real(s(:)))+10]);
    title('Real part');
    subplot(412);plot(te(:).',imag(permute(s,[2 1])),'k^-');hold on;ylim([min(imag(s(:)))-10,max(imag(s(:)))+10]);
    title('Imaginary part');
    subplot(413);plot(te(:).',abs(permute(s,[2 1])),'k^-');hold on;ylim([0,max(abs(s(:)))+10]);
    title('Magnitude');
    subplot(414);plot(te(:).',angle(permute(s,[2 1])),'k^-');hold on;ylim([-4 4]);
    title('Phase');
    subplot(411);plot(te(:).',real(permute(sHat,[2 1])),'x-');plot(te(:).',real(permute(sHat,[2 1]))-real(permute(s,[2 1])),'ro-.');
    hold off;
    subplot(412);plot(te(:).',imag(permute(sHat,[2 1])),'x-');plot(te(:).',imag(permute(sHat,[2 1]))-imag(permute(s,[2 1])),'ro-.');
    hold off;
    subplot(413);plot(te(:).',abs(permute(sHat,[2 1])),'x-');plot(te(:).',abs(abs(permute(sHat,[2 1]))-abs(permute(s,[2 1]))),'ro-.');
    hold off;
    text(te(1)*1.1,max(abs(s(:))*0.9),sprintf('resnorm=%f',sum(err(:).^2)));
    text(te(1)*1.1,max(abs(s(:))*0.8),sprintf('Amy=%f,Aax=%f,Aex=%f,t2*my=%f,t2*ax=%f,t2*ex=%f,fmy=%f,fax=%f,fex=%f,pini=%f',...
        Amy,Aax,Aex,t2smy,t2sax,t2sex,fmybg,faxbg,fexbg,pini));
    subplot(414);plot(te(:).',angle(permute(sHat,[2 1])),'x-');plot(te(:).',angle(permute(s.*conj(sHat),[2 1])),'ro-.');
    hold off;
    drawnow;
end

end

%% check and set default
function [algoPara2,imgPara2,isValid]=CheckAndSetPara(algoPara,imgPara)

imgPara2 = imgPara;
algoPara2 = algoPara;
isValid = true;

% check if the number of echo times matches with the data
if length(imgPara.te) ~= size(imgPara.img,4)
    isValid = false;
end

% check debug
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

% check signal mask
try
    imgPara2.mask = imgPara.mask;
catch
    imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
end

end