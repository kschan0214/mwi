load_module_mwi;
clear;
%%
% varying parameters
echoSpacing = [1:0.2:5]*1e-3;
nTE = [12:4:128];
% TR_r      = [50:10:80]*1e-3;
% nSpacing = length(echoSpacing);
te1 = 2.15e-3;
% te_r = struct();
% 
% for kte = 1:length(nTE)
%     for kSpacing = 1:nSpacing
%         cte = te1:echoSpacing(kSpacing):(nTE(kte)-1)*echoSpacing(kSpacing);
%         te_r(kte,kSpacing).te = cte;
%     end
% end
% varying parameters
TR_r      = [50:10:500]*1e-3;
nSpacing = length(echoSpacing);
te_r = struct();
for ktr = 1:length(TR_r)
    ctr = TR_r(ktr);
    for kSpacing = 1:nSpacing
        for kt = 1:200
            curr_te = te1 + (kt-1)*echoSpacing(kSpacing);
            if curr_te < ctr
                te_r(ktr,kSpacing).te(kt) = curr_te;
            end
        end
    end
end
%%
mwf = 12/100;
icvf = 60/90;
S0 = 100;

t2sMyelin = 10e-3;   %s
t2sIW= 64e-3;    %s
t2sEW = 48e-3;    %s

freqMyelin = [6.70928162780528];
freqIW = [-1.56906557957569];
freqEW = 0;	%Hz

r2sA = 1/t2sMyelin;
r2sB = 1/t2sIW;
r2sC = 1/t2sEW;

wA = 2*pi*freqMyelin;
wB = 2*pi*freqIW;
wC = 2*pi*freqEW;

phi = 0;
w_mb = wA;
w_ib = wB;
w_eb = wC;
r2s_mw = r2sA;
r2s_iw = r2sB;
r2s_ew = r2sC;

% fim_all = zeros(10,10,length(nTE),nSpacing);
% crlb_all = zeros(10,length(nTE),nSpacing);
fim_all3 = zeros(7,7,length(TR_r),nSpacing);
crlb_all3 = zeros(7,length(TR_r),nSpacing);
sigma = 0.5;
for kte = 1:length(TR_r)
    for kSpacing = 1:nSpacing
        t = te_r(kte,kSpacing).te;
        [fim, crlb] = mwi_crlb_3cc_mwf_GivenicvfHCFM(t,sigma,S0,mwf,icvf,r2s_mw,r2s_iw,r2s_ew,w_mb,w_ib,w_eb,phi);
        fim_all3(:,:,kte,kSpacing) = fim;
        crlb_all3(:,kte,kSpacing) = crlb;
    end
end
