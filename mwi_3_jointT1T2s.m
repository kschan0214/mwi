%% fitRes = mwi_3_jointT1T2s(algoPara,imgPara)
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
% Date created: 4 january 2018
% Date last modified:
%
%
function fitRes = mwi_3_jointT1T2s(algoPara,imgPara)

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

% capture all parameters
numMagn    = algoPara.numMagn;
maxIter    = algoPara.maxIter;
isParallel = algoPara.isParallel;

te    = imgPara.te;
tr    = imgPara.tr;
fa    = imgPara.fa;
data  = imgPara.img;
mask  = imgPara.mask;
b0map = imgPara.b0map;
b1map = imgPara.b1map;

[nx,ny,nz,~,~] = size(data);

options = optimoptions(@lsqnonlin,'Display','off','MaxIter',maxIter);

estimates = zeros(nx,ny,nz,13);
resnorm   = zeros(nx,ny,nz);
for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            if mask(kx,ky,kz)>0
                s = permute(data(kx,ky,kz,:,:),[5 4 1 2 3]);
                b0 = b0map(kx,ky,kz);
                b1 = b1map(kx,ky,kz);
                if kx==76 && ky==22
                [estimates(kx,ky,kz,:),resnorm(kx,ky,kz)] = FitModel(s,fa,te,tr,b0,b1,numMagn,options);
                end
            end
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b0,b1,numMagn,options)
% set initial guesses
[t10,m0] = DESPOT1(abs(s(:,1)),fa,tr);
[S1max,ind] = max(abs(s(:,1)));
Amy0   = 0.1*m0;            Amylb   = 0;        Amyub   = 2*m0;
Aax0   = 0.6*m0;            Aaxlb   = 0;        Aaxub   = 2*m0;
Aex0   = 0.3*m0;            Aexlb   = 0;        Aexub   = 2*m0;
t2smy0 = 10e-3;             t2smylb = 3e-3;     t2smyub = 25e-3;
t2sax0 = 64e-3;             t2saxlb = 25e-3;    t2saxub = 150e-3;
t2sex0 = 48e-3;             t2sexlb = 25e-3;    t2sexub = 150e-3;
t1my0  = 300e-3;            t1mylb  = 150e-3;   t1myub  = 500e-3;
% t1ax0  = 1000e-3;           t1axlb  = 500e-3;   t1axub  = 1600e-3;
% t1ex0  = 1200e-3;           t1exlb  = 500e-3;   t1exub  = 1600e-3;
t1ax0  = t10;               t1axlb  = 500e-3;   t1axub  = 1800e-3;
t1ex0  = t10+0.2;           t1exlb  = 500e-3;   t1exub  = 1800e-3;
fmy0   = b0;                fmylb   = b0-75;    fmyub   = b0+75;
fax0   = b0;                faxlb   = b0-25;    faxub   = b0+25;
fex0   = b0;                fexlb   = b0-25;    fexub   = b0+25;
pini0  = angle(s(1,ind));   pinilb  = -pi;      piniub  = pi;

x0 = [Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,t1my0,t1ax0,t1ex0,fmy0,fax0,fex0,pini0];
lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,t1mylb,t1axlb,t1exlb,fmylb,faxlb,fexlb,pinilb];
ub = [Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,t1myub,t1axub,t1exub,fmyub,faxub,fexub,piniub];

[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,numMagn),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,numMagn)
Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5); t2sex=x(6);
t1my=x(7);  t1ax=x(8);  t1ex=x(9);
fmy=x(10);  fax=x(11);  fex=x(12);
pini=x(13);

sHat = mwi_model_3cc_jointT1T2s(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1ax,t1ex,fmy,fax,fex,pini,b1);

err = computeFiter(s,sHat,numMagn);

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

% check number of phase-corrupted echoes
try
    algoPara2.numMagn = algoPara.numMagn;
catch
    algoPara2.numMagn = 0;
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
% check b0 map
try
    imgPara2.b0map = imgPara.b0map;
catch
    imgPara2.b0map = zeros(size(imgPara.img,1),size(imgPara.img,2),size(imgPara.img,3));
end
% check b1 map
try
    imgPara2.b1map = imgPara.b1map;
catch
    imgPara2.b1map = ones(size(imgPara.img,1),size(imgPara.img,2),size(imgPara.img,3));
end

end