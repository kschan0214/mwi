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
function fitRes = mwi_3cm_jointT1T2s_DEBUG(algoPara,imgPara)
global DEBUG
DEBUG=1;

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

% capture all parameters
numMagn    = length(imgPara.te);
maxIter    = algoPara.maxIter;
% isParallel = algoPara.isParallel;

te    = imgPara.te;
tr    = imgPara.tr;
fa    = imgPara.fa;
data  = imgPara.img;
mask  = imgPara.mask;
b1map = imgPara.b1map;

[ny,nx,nz,~,~] = size(data);

if DEBUG
%     fdss = [1e-2,1e-2,1e-3,1e-7,1e-7,1e-7,1e-6,1e-6,1e-6,1e-5,1e-5];
    fdss = [1e-4,1e-4,1e-4,1e-7,1e-7,1e-7,1e-6,1e-6,1e-5,1e-5];
    options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'FiniteDifferenceStepSize',fdss,'MaxFunctionEvaluations',200*11);
else
    options = optimoptions(@lsqnonlin,'Display','off','MaxIter',maxIter);
end

estimates = zeros(ny,nx,nz,10);
resnorm   = zeros(ny,nx,nz);
for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            if mask(ky,kx,kz)>0
                % 1st dim: T1w; 2nd dim: T2*w
                s = permute(data(ky,kx,kz,:,:),[5 4 1 2 3]);
                b1 = b1map(ky,kx,kz);
                if ky==38 && kx==63 && DEBUG
                [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,fa,te,tr,b1,numMagn,options);
                end
            end
        end
%         if DEBUG && ky==round(ny/2)
%             figure;imshow(estimates(:,:,kz,4),[0 0.02]);
%             drawnow;
%         end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,fa,te,tr,b1,numMagn,options)
global DEBUG
% set initial guesses
t2s0 = zeros(1,length(fa));
s0 = zeros(1,length(fa));
for kt=1:length(te)-1
    for kfa=1:length(fa)
        [~,t2s0(kfa,kt),s0(kfa,kt)] = R2starmapping_trapezoidal_voxel(s(kfa,kt:kt+1),te(kt:kt+1));
    end
end
t100 = zeros(1,length(te));
m00 = zeros(1,length(te));
for kt=1:length(te)
    [t100(kt),m00(kt)] = DESPOT1(squeeze(abs(s(:,kt))),fa,tr);
end

% [t10,m0] = DESPOT1(abs(s0(:,1)),fa,tr);

[t10,m0] = DESPOT1(abs(s(:,1)),fa,tr);

% [t10,m0] = DESPOT1(abs(s(:,1)),fa,tr);
if t10>2500e-3 || t10<500e-3
    t10=1000e-3;
end
if m0<0
    m0=0;
end
% [S1max,ind] = max(abs(s(:,1)));
Amy0   = 0.1*m0;            Amylb   = 0;        Amyub   = 0.2*m0;
Aax0   = 0.6*m0;            Aaxlb   = 0;        Aaxub   = 1*m0;
Aex0   = 0.3*m0;            Aexlb   = 0;        Aexub   = 1*m0;
% % Aax0   = 0.7*m0;            Aaxlb   = 0.5*m0;        Aaxub   = 1*m0;
% % Aex0   = 0.3*m0;            Aexlb   = 0.2*m0;        Aexub   = 0.5*m0;
% t2smy0 = 10e-3;             t2smylb = 3e-3;     t2smyub = 25e-3;
% t2sax0 = t2s0(1);             t2saxlb = 25e-3;    t2saxub = 150e-3;
% t2sex0 = t2s0(1)-10e-3;             t2sexlb = 25e-3;    t2sexub = 150e-3;
% t1my0  = 100e-3;            t1mylb  = 150e-3;   t1myub  = 500e-3;
% % t1ax0  = 1000e-3;           t1axlb  = 500e-3;   t1axub  = 1600e-3;
% % t1ex0  = 1200e-3;           t1exlb  = 500e-3;   t1exub  = 1600e-3;
% t1ax0  = t10;               t1axlb  = 500e-3;   t1axub  = 1500e-3;
% t1ex0  = t10-0.2;           t1exlb  = 500e-3;   t1exub  = 1500e-3;
% fmy0   = 11;                fmylb   = 11-75;    fmyub   = 11+75;
% fax0   = 0;                 faxlb   = -25;      faxub   = +25;
% Amy0   = 0.1*m00(end);            Amylb   = 0;        Amyub   = 0.2*m00(end);
% Aax0   = 0.6*m00(end);            Aaxlb   = 0;        Aaxub   = 1*m00(end);
% Aex0   = 0.3*m00(end);            Aexlb   = 0;        Aexub   = 1*m00(end);
% Aax0   = 0.7*m0;            Aaxlb   = 0.5*m0;        Aaxub   = 1*m0;
% Aex0   = 0.3*m0;            Aexlb   = 0.2*m0;        Aexub   = 0.5*m0;
t2smy0 = 10e-3;             t2smylb = 3e-3;     t2smyub = 25e-3;
t2sax0 = 64e-3;      t2saxlb = 25e-3;    t2saxub = 150e-3;
t2sex0 = 48e-3;             t2sexlb = 25e-3;    t2sexub = 150e-3;
% t2sl0 = max(t2s0(:));             t2sllb = 25e-3;    t2slub = 150e-3;
% t1ax0  = 1000e-3;           t1axlb  = 500e-3;   t1axub  = 1600e-3;
% t1ex0  = 1200e-3;           t1exlb  = 500e-3;   t1exub  = 1600e-3;
t1my0  = 200e-3;            t1mylb  = 50e-3;   t1myub  = 500e-3;
t1l0  = mean(t100);               t1llb  = 500e-3;   t1lub  = 2000e-3;
fmy0   = 5;                fmylb   = -10;    fmyub   = 5+75;
fax0   = 0;                 faxlb   = -25;      faxub   = 10;
% b10 = 1; b1lb=0; b1ub=1.5;
% x0 = [Amy0,Aax0,Aex0,t2smy0,t2sl0,t1my0,t1l0,fmy0,fax0,b10];
% lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2sllb,t1mylb,t1llb,fmylb,faxlb,b1lb];
% ub = [Amyub,Aaxub,Aexub,t2smyub,t2slub,t1myub,t1lub,fmyub,faxub,b1ub];
x0 = [Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,t1my0,t1l0,fmy0,fax0];
lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,t1mylb,t1llb,fmylb,faxlb];
ub = [Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,t1myub,t1lub,fmyub,faxub];

% x0 = [Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,t1my0,t1ax0,t1ex0,fmy0,fax0];
% lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,t1mylb,t1axlb,t1exlb,fmylb,faxlb];
% ub = [Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,t1myub,t1axub,t1exub,fmyub,faxub];

[x,res] = lsqnonlin(@(y)CostFunc(y,s,fa,te,tr,b1,numMagn),x0,lb,ub,options);

if DEBUG
    s_avg = abs(s(1,:));
    [~,t2s0_debug,s0_debug] = R2starmapping_trapezoidal_voxel(s_avg,te);
    s_simulate = Signal_mGRE(s0_debug,t2s0_debug,te);
    t2s_norm = sum(sum(bsxfun(@minus,abs(s),abs(s_simulate)).^2));
    fprintf('T2*-only residual norm = %f\n',t2s_norm);

    [t10_debug,m0_debug] = DESPOT1(abs(s(:,1)),fa,tr);
    S_simulated = Signal_GRE_T1wMono(m0_debug, fa, t10_debug, tr);
    t1_norm = sum((abs(s(:,1))-S_simulated).^2);

    s_simulated = Signal_GRE_mcJointT1T2star(m0,t10,fa,tr,mean(t2s0),te,0);
    normAll = sum((abs(s(:))-s_simulated(:)).^2);
    figure(98);
%     subplot(121);
    plot(abs(permute(s,[2 1])),'k');hold on;
    plot(abs(permute(s_simulated,[2 1])),'x');
%     hold off;
%     title('1-comp');
    fprintf('1-compartment residual norm = %f\n',normAll);
end
if DEBUG
%     s_simulated2 = mwi_model_3cm_jointT1T2s(fa,te,tr,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11));
    s_simulated2 = mwi_model_3cm_jointT1T2s(fa,te,tr,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(8),x(9),x(10));
    figure(98);
%     subplot(121);
    plot(abs(permute(s,[2 1])),'k');hold on;
    plot(abs(permute(s_simulated2,[2 1])),'^');hold off;
    title('3-comp vs 1-comp');
    fprintf('3-compartment residual norm = %f\n',sum((abs(s(:))-s_simulated2(:)).^2));
end


end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,fa,te,tr,b1,numMagn)
global DEBUG
% Amy=x(1);   Aax=x(2);   Aex=x(3);
% t2smy=x(4); t2sax=x(5); t2sex=x(6);
% t1my=x(7);  t1ax=x(8);  t1ex=x(9);
% fmy=x(10);  fax=x(11);

Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5);  t2sex=x(6);
t1my=x(7);  t1l=x(8);
fmy=x(9);  fax=x(10);

% b1=x(10);
% sHat = mwi_model_3cm_jointT1T2s(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1ax,t1ex,fmy,fax,b1);

% sHat = mwi_model_3cm_jointT1T2s_2(fa,te,tr,Amy,Aax,Aex,t2smy,t2sl,t1my,t1l,fmy,fax,b1);

sHat = mwi_model_3cm_jointT1T2s_2(fa,te,tr,Amy,Aax,Aex,t2smy,t2sax,t2sex,t1my,t1l,fmy,fax,b1);

err = computeFiter(s,sHat,numMagn);

if DEBUG
%     w = abs(s(:))/norm(abs(s(:)));
%         w = ones(size(s));
%         w(1,4,:) = 100;
%         w = bsxfun(@times,abs(s),1./vecnorm(abs(s),2,2));
%         w = w(:)./norm(w(:));
%     err = err.*w(:);
    figure(99);plot(abs(permute(s,[2 1])),'k');hold on;ylim([min(abs(s(:)))-10,max(abs(s(:)))+10]);
    plot(abs(permute(sHat,[2 1])),'x-');hold off;text(size(s,1)/2,max(abs(s(:))),sprintf('resnorm = %f',sum(err.^2)));
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