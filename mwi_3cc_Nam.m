%% fitRes = mwi_3cc_Nam(algoPara,imgPara)
%
% Input
% --------------
% algoPara.maxIter : maximum iteration allows (default 500)
% algoPara.isROI   : boolean ROI analysis (default false)
% algoPara.DEBUG   : debug mode (default false)
% algoPara.fcnTol  : function tolerance (default: 5 orders small than max. signal)
% imgPara.img      : 4D image data, time in 4th dimension
% imgPara.mask     : signal mask
% imgPara.te       : echo times
% imgPara.fieldmap : background field
%
% Output
% --------------
% fitRes.estimates : fitting estimates (Ampl_n,t2s_n,freq_n)
% fitres.resnorm   : L2 norm of fitting residual
%
% Description: Myelin water mapping by fitting complex mode(c) with
% complex-valued data(c) (ref. Model 3 in Nam et al. 2015 NeuroImage)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 19 January 2018
% Date last modified: 27 February 2018
%
%
function fitRes = mwi_3cc_Nam(algoPara,imgPara)

% check validity of the algorithm parameters and image parameters
[algoPara,imgPara,isValid]=CheckAndSetPara(algoPara,imgPara);
if ~isValid
    fitRes = [];
    disp('Invalid parameters');
    return
end

DEBUG   = algoPara.DEBUG;
verbose = algoPara.verbose;

% capture all parameters
numMagn    = 0;
maxIter    = algoPara.maxIter;
fcnTol     = algoPara.fcnTol;
isWeighted = algoPara.isWeighted;
% isROI      = algoPara.isROI; 
% isParallel = algoPara.isParallel;

te    = imgPara.te;
data  = imgPara.img;
mask  = imgPara.mask;
fm    = imgPara.fieldmap;

% for complex fitting we have doubled the elements in the cost function
fcnTol = fcnTol./(2*numel(te));

% display fitting message
if verbose
    disp('The following fitting parameters are used:');
    fprintf('# max. iteration=%i\n',maxIter);
    fprintf('Function tolerance=%e\n',fcnTol);
    if isWeighted
        disp('Cost function is weighted by echo intensity');
    else
        disp('No weight on cost function');
    end
end

[ny,nx,nz,~,~] = size(data);

% get the order of magnitude of signal
orderS = floor(log10(max(abs(data(:)))));

options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*15,...
    'StepTolerance',1e-6,'FunctionTolerance',fcnTol);

% fdss = [10^(orderS-6),10^(orderS-6),10^(orderS-6),1e-8,1e-8,1e-8,1e-8,1e-8,1e-8],1e-8]];
% options.FiniteDifferenceStepSize = fdss;
if ~DEBUG
    options.Display = 'off';
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
                [estimates(ky,kx,kz,:),resnorm(ky,kx,kz)] = FitModel(s,te,db0,numMagn,isWeighted,options,DEBUG);
            end
        end
    end
end

fitRes.estimates = estimates;
fitRes.resnorm   = resnorm;

end

%% Setup lsqnonlin and fit with the default model
function [x,res] = FitModel(s,te,db0,numMagn,isWeighted,options,DEBUG)
if DEBUG
    % if DEBUG then create an array to store resnorm of all iterations
    global DEBUG_resnormAll
    DEBUG_resnormAll=[];
end

% [~,~,s0]=R2starmapping_trapezoidal_voxel(abs(s),te);
Amy0   = 0.1*abs(s(1));            Amylb   = 0;        Amyub   = 1*abs(s(1));
Aax0   = 0.6*abs(s(1));            Aaxlb   = 0;        Aaxub   = 1*abs(s(1));
Aex0   = 0.3*abs(s(1));            Aexlb   = 0;        Aexub   = 1*abs(s(1));
t2smy0 = 10e-3;                    t2smylb = 3e-3;     t2smyub = 25e-3;
t2sax0 = 64e-3;               t2saxlb = 25e-3;    t2saxub = 350e-3;
t2sex0 = 48e-3;               t2sexlb = 25e-3;    t2sexub = 350e-3;
fmy0   = db0;                 fmylb   = db0-75;    fmyub   = db0+75;
fax0   = db0;                 faxlb   = db0-25;      faxub   = db0+25;
fex0   = db0;                 fexlb   = db0-25;      fexub   = db0+25;
pini0   = 2*pi*db0*te(1)-angle(s(1));        pinilb = -pi;         piniub=pi;


x0 = double([Amy0,Aax0,Aex0,t2smy0,t2sax0,t2sex0,fmy0,fax0,fex0,pini0]);
lb = double([Amylb,Aaxlb,Aexlb,t2smylb,t2saxlb,t2sexlb,fmylb,faxlb,fexlb,pinilb]);
ub = double([Amyub,Aaxub,Aexub,t2smyub,t2saxub,t2sexub,fmyub,faxub,fexub,piniub]);

[x,res] = lsqnonlin(@(y)CostFunc(y,s,te,numMagn,isWeighted,DEBUG),x0,lb,ub,options);

end

%% compute the cost function of the optimisation problem
function err = CostFunc(x,s,te,numMagn,isWeighted,DEBUG)
% obatin fitting parameters
Amy=x(1);   Aax=x(2);   Aex=x(3);
t2smy=x(4); t2sax=x(5); t2sex=x(6);
fmybg=x(7);  faxbg=x(8);    fexbg=x(9);
pini=x(10);

% simulate signal based on parameter input
sHat = mwi_model_3cc_nam2015(te,Amy,Aax,Aex,t2smy,t2sax,t2sex,fmybg,faxbg,fexbg,pini);

% compute fitting residual
if isWeighted
    % weighted the cost function by echo intensity, as suggested in Nam's paper
    w = abs(s(:))/norm(abs(s(:)));
    err = computeFiter(s,sHat,numMagn,w);
%     w = abs(s(:));
%     err = err.*w(:);
else
    err = computeFiter(s,sHat,numMagn);
end

% cost function is normalised with the norm of signal in order to provide
% sort of consistence with fixed function tolerance
err = err ./ norm(abs(s(:)));

if DEBUG
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
% check parallel computing 
try
    algoPara2.isParallel = algoPara.isParallel;
catch
    algoPara2.isParallel = false;
end
% check function tolerance
try
    algoPara2.fcnTol = algoPara.fcnTol;
catch
    % get the order of magnitude of signal
%     orderS = floor(log10(max(abs(imgPara2.img(:)))));
%     algoPara2.fcnTol = 10^(orderS-5);
    algoPara2.fcnTol = 1e-5;
end
% check weighted sum of cost function
try
    algoPara2.isWeighted = algoPara.isWeighted;
catch
    algoPara2.isWeighted = true;
end

% check ROI analysis
try
    algoPara2.isROI = algoPara.isROI;
catch
    algoPara2.isROI = false;
end

% check signal mask
try
    imgPara2.mask = imgPara.mask;
catch
    imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
end

end