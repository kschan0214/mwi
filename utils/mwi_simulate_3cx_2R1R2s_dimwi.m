%% function output = function_name(input)
%
% Usage:
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
% Date created: 
% Date last modified:
%
%

function img_hat = mwi_simulate_3cx_2R1R2s_dimwi(algoPara,imgPara,fitres)

if isfield(fitres,'mask')
    mask = fitres.mask;
elseif isfield(imgPara,'mask')
    mask = imgPara.mask;
end

b0dir = double(imgPara.b0dir);

freq_ew = 0;

DIMWI = algoPara.DIMWI;
% store fixed parameters
fixParam.b0     = double(imgPara.b0);
fixParam.rho_mw = double(imgPara.rho_mw);
fixParam.E      = double(imgPara.E);
fixParam.x_i    = double(imgPara.x_i);
fixParam.x_a    = double(imgPara.x_a);
DIMWI.fixParam = fixParam;


epgx.npulse       = algoPara.npulse;
epgx.rfphase      = algoPara.rfphase;
epgx.isExchange   = algoPara.isExchange;
epgx.isT1mw       = algoPara.isT1mw;
epgx.isEPG        = algoPara.isEPG;
epgx.T1mw         = algoPara.T1mw;
epgx.rho_mw       = fixParam.rho_mw;

% check if fibre orientation is needed
if DIMWI.isFreqMW || DIMWI.isFreqIW || DIMWI.isR2sEW
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
else
    ff      = ones(size(mask));
    theta   = zeros(size(mask));
end

tr = imgPara.tr;
fa = imgPara.fa;
te = imgPara.te;

img_hat = zeros(size(imgPara.img));
for kz = 1:size(fitres.estimates,3)
    for ky = 1:size(fitres.estimates,2)
        for kx = 1:size(fitres.estimates,1)
            if mask(kx,ky,kz) == 1
                
                b1      = imgPara.b1map(kx,ky,kz);
                
                Amw     = fitres.S0_MW(kx,ky,kz);
                if DIMWI.isVic
                    DIMWI.icvf  = imgPara.icvf(kx,ky,kz);
                    Aiw = fitres.S0_IEW * DIMWI.icvf;
                    Aew = fitres.S0_IEW * (1-DIMWI.icvf);
                else
                    Aiw = fitres.S0_IW(kx,ky,kz);
                    Aew = fitres.S0_EW(kx,ky,kz);
                end
                
                t2smw   = 1/fitres.R2s_MW(kx,ky,kz);
                t2siw   = 1/fitres.R2s_IW(kx,ky,kz);
                if DIMWI.isR2sEW 
                    t2sew = 1/fitres.R2s_IW(kx,ky,kz);
                else
                    t2sew   = 1/fitres.R2s_EW(kx,ky,kz);
                end
                if DIMWI.isFreqMW  
                    freq_mw = 0;
                else
                    freq_mw = fitres.Freq_MW(kx,ky,kz);
                end
                if DIMWI.isFreqIW
                    freq_iw = 0;
                else
                    freq_iw = fitres.Freq_IW(kx,ky,kz);
                end
                
                fbg     = fitres.Freq_BKG(kx,ky,kz,:);
                pini    = fitres.pini(kx,ky,kz,:);
                t1l     = fitres.T1_IEW(kx,ky,kz);
                
                if algoPara.isT1mw
                    t1s = fitres.T1_M(kx,ky,kz);
                else
                    t1s = algoPara.T1mw;
                end
                
                if algoPara.isExchange == 1
                    kls = fitres.kiewm;
                else
                    kls = 0;
                end
                
                if epgx.isExchange && epgx.isEPG
                    % precompute EPG-X's transition matrix here for speed
                    phiCycle = RF_phase_cycle(epgx.npulse,epgx.rfphase);
                    for kfa=1:length(fa)
                        T3D_all{kfa} = PrecomputeT(phiCycle,d2r(fa(kfa)*b1));
                    end
                    epgx.T3D_all = T3D_all;
                end
                
                
                if DIMWI.isFreqMW || DIMWI.isFreqIW || DIMWI.isR2sEW
                    DIMWI.theta = theta(kx,ky,kz,:);
                    DIMWI.ff    = ff(kx,ky,kz);
                else
                    DIMWI.theta = 0;
                    DIMWI.ff    = 1;
                end
                
                sHat = zeros([length(te), length(fa), length(DIMWI.theta)]);
                DIMWI_curr = DIMWI;
                for kfo = 1:size(DIMWI.theta,4)
    
                    DIMWI_curr.theta = DIMWI.theta(kfo);

                    sHat(:,:,kfo) = mwi_model_2T13T2scc_dimwi(te,tr,fa,b1,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1s,t1l,kls,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI_curr,epgx);
                end
                img_hat(kx,ky,kz,:,:) = sum(bsxfun(@times,sHat,permute(DIMWI.ff(:),[2 3 1])),3);
            end
        end
    end
end

end