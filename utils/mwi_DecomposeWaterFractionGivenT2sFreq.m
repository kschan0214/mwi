function A = mwi_DecomposeWaterFractionGivenT2sFreq(s,te,t2s,freq,totalField,pini)
if nargin<5
    totalField = 0;
end
if nargin<6
    pini = 0;
end

ncomp = length(t2s);

relaxation = zeros(1,ncomp);
for kcomp=1:ncomp
    relaxation(kcomp) = -1/t2s(kcomp)+1i*2*pi*freq(kcomp);
end
relaxation = exp(te(:) * relaxation);
external_field = diag(exp(1i*2*pi*totalField*te) * exp(-1i*pini));

% A = pinv(external_field*relaxation) * s(:);
A = (external_field*relaxation) \ s(:);

end