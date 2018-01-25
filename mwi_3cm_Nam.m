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
data  = imgPara.img;
mask  = imgPara.mask;

[ny,nx,nz,~,~] = size(data);

if DEBUG
%     fdss = [1e-2,1e-2,1e-3,1e-7,1e-7,1e-7,1e-6,1e-6,1e-6,1e-5,1e-5];
    fdss = [1e-4,1e-4,1e-4,1e-7,1e-7,1e-7,1e-5,1e-5];
    options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'FiniteDifferenceStepSize',fdss,'MaxFunctionEvaluations',200*15,...
        'StepTolerance',1e-5,'FunctionTolerance',1e-5);
else
    options = optimoptions(@lsqnonlin,'Display','off','MaxIter',maxIter);
end

estimates = zeros(ny,nx,nz,8);
resnorm   = zeros(ny,nx,nz);
for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            if mask(ky,kx,kz)>0
                % 1st dim: T1w; 2nd dim: T2*w
                s = permute(data(ky,kx,kz,:),[5 4 1 2 3]);
                if ky==47 && kx==100 && DEBUG
                [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,numMagn,options);
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
function [x,res] = FitModel(s,te,numMagn,options)
global DEBUG
% set initial guesses
t2s0 = zeros(1,length(te)-1);
s0 = zeros(1,length(te)-1);
for kt=1:length(te)-1
    [~,t2s0(kt),s0(kt)] = R2starmapping_trapezoidal_voxel(s(kt:kt+1),te(kt:kt+1));
end

% [S1max,ind] = max(abs(s(:,1)));

% Amy0   = 0.1*max(s0);            Amylb   = 0;        Amyub   = 0.2*max(s0);
% Aax0   = 0.6*max(s0);            Aaxlb   = 0;        Aaxub   = 1*max(s0);
% Aex0   = 0.3*max(s0);            Aexlb   = 0;        Aexub   = 1*max(s0);
Amy0   = 0.1*abs(s(1));            Amylb   = 0;        Amyub   = 0.3*abs(s(1));
Aax0   = 0.6*abs(s(1));            Aaxlb   = 0;        Aaxub   = 1*abs(s(1));
Aex0   = 0.3*abs(s(1));            Aexlb   = 0;        Aexub   = 1*abs(s(1));
t2smy0 = 5e-3;             t2smylb = 1e-3;     t2smyub = 20e-3;
t2sax0 = max(t2s0(:));             t2saxlb = 20e-3;    t2saxub = 150e-3;
t2sex0 = max(t2s0(:))-10e-3;             t2sexlb = 20e-3;    t2sexub = 150e-3;
fmy0   = 5;                fmylb   = 5-75;    fmyub   = 5+75;
fax0   = 0;                 faxlb   = -25;      faxub   = +25;


x0 = [Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,fmy0,fax0];
lb = [Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,fmylb,faxlb];
ub = [Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,fmyub,faxub];


[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,numMagn),x0,lb,ub,options);


% if DEBUG
%     s_avg = abs(s(1,:));
%     [~,t2s0_debug,s0_debug] = R2starmapping_trapezoidal_voxel(s_avg,te);
%     s_simulate = Signal_mGRE(s0_debug,t2s0_debug,te);
%     t2s_norm = sum(sum(bsxfun(@minus,abs(s),abs(s_simulate)).^2));
%     fprintf('T2*-only residual norm = %f\n',t2s_norm);
% 
%     [t10_debug,m0_debug] = DESPOT1(abs(s(:,1)),fa,tr);
%     S_simulated = Signal_GRE_T1wMono(m0_debug, fa, t10_debug, tr);
%     t1_norm = sum((abs(s(:,1))-S_simulated).^2);
% 
%     s_simulated = Signal_GRE_mcJointT1T2star(m0,t10,fa,tr,mean(t2s0),te,0);
%     normAll = sum((abs(s(:))-s_simulated(:)).^2);
%     figure(98);
% %     subplot(121);
%     plot(abs(permute(s,[2 1])),'k');hold on;
%     plot(abs(permute(s_simulated,[2 1])),'x');
% %     hold off;
% %     title('1-comp');
%     fprintf('1-compartment residual norm = %f\n',normAll);
% end
% if DEBUG
% %     s_simulated2 = mwi_model_3cm_jointT1T2s(fa,te,tr,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11));
%     s_simulated2 = mwi_model_3cm_jointT1T2s(fa,te,tr,x(1),x(2),x(3),x(4),x(5),x(5),x(6),x(7),x(7),x(8),x(9));
%     figure(98);
% %     subplot(121);
%     plot(abs(permute(s,[2 1])),'k');hold on;
%     plot(abs(permute(s_simulated2,[2 1])),'^');hold off;
%     title('3-comp vs 1-comp');
%     fprintf('3-compartment residual norm = %f\n',sum((abs(s(:))-s_simulated2(:)).^2));
% end


end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,numMagn)
global DEBUG
Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5); t2sex=x(6);
fmy=x(7);  fax=x(8);

sHat = mwi_model_3cm_vGelderen2012(te,Amy,Aax,Aex,t2smy,t2sax,t2sex,fmy,fax);

err = computeFiter(s,sHat,numMagn);

if DEBUG
%     w = abs(s(:))/norm(abs(s(:)));
%         w = ones(size(s));
%         w(1,:) = 2;
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

end