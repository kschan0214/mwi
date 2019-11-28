
function [freq_myelin, freq_axon] = hcfm_freq_givenDIMWI(fitres,imgPara)

ff    = double(imgPara.ff); % fibre fraction
try 
    fo    = double(imgPara.fo); % fibre orientation
    theta = zeros(size(ff));
    for kfo = 1:size(fo,5)
        theta(:,:,:,kfo) = AngleBetweenV1MapAndB0(fo(:,:,:,:,kfo),imgPara.b0dir);
    end
catch 
    theta = double(imgPara.theta); % theta map
end
% normalise fibre fraction
ff = bsxfun(@rdivide,ff,sum(ff,4));
ff(isnan(ff)) = 0;

Amw = fitres.S0_MW;
Aiw = fitres.S0_IEW .* imgPara.icvf;

g = sqrt(abs(Aiw)./(abs(Aiw)+abs(Amw)/imgPara.rho_mw));

freq_myelin = zeros(size(ff));
freq_axon   = zeros(size(ff));
for kfo = 1:size(theta,4)

    [freq_myelin(:,:,:,kfo), freq_axon(:,:,:,kfo)] = hcfm_freq(imgPara.x_i,imgPara.x_a,g,theta(:,:,:,kfo),imgPara.E,imgPara.b0);

end
freq_myelin = sum(freq_myelin.*ff,4);
freq_axon   = sum(freq_axon.*ff,4);

freq_myelin(isnan(freq_myelin)) = 0;
freq_axon(isnan(freq_axon)) = 0;
freq_myelin(isinf(freq_myelin)) = 0;
freq_axon(isinf(freq_axon)) = 0;

end