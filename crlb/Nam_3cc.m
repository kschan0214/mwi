load_module_mwi;
clear;
%%
% varying parameters
echoSpacing = [1:0.1:5]*1e-3;
nTE = [10:1:64];
% TR_r      = [50:10:80]*1e-3;
nSpacing = length(echoSpacing);
te1 = 2.15e-3;
te_r = struct();

for kte = 1:length(nTE)
    for kSpacing = 1:nSpacing
        cte = te1:echoSpacing(kSpacing):te1+(nTE(kte)-1)*echoSpacing(kSpacing);
        te_r(kte,kSpacing).te = cte;
    end
end
%% varying parameters
% TR_r      = [50:10:400]*1e-3;
% nSpacing = length(echoSpacing);
% te_r = struct();
% for ktr = 1:length(TR_r)
%     ctr = TR_r(ktr);
%     for kSpacing = 1:nSpacing
%         for kt = 1:500
%             curr_te = te1 + (kt-1)*echoSpacing(kSpacing);
%             if curr_te < ctr
%                 te_r(ktr,kSpacing).te(kt) = curr_te;
%             end
%         end
%     end
% end
%%
mwf     = 12/100;
icvf    = 60/90;
S0      = 100;

t2sMyelin   = 10e-3;   %s
t2sIW       = 64e-3;    %s
t2sEW       = 48e-3;    %s

g       = 0.8;
theta   = 90;
x_i     = -0.1;
x_a     = -0.1;
E       = 0.02;
b0      = 3;
[freqMyelin,freqIW] = HCFM_FreqShift(theta,x_i,x_a,g,E,b0);
% freqMyelin = [6.70928162780528];
% freqIW = [-1.56906557957569];
freqEW = 0;	%Hz

phi     = 0;
w_mb    = 2*pi*freqMyelin;
w_ib    = 2*pi*freqIW;
w_eb    = 2*pi*freqEW;
r2s_mw  = 1/t2sMyelin;
r2s_iw  = 1/t2sIW;
r2s_ew  = 1/t2sEW;

fim_all = zeros(10,10,length(nTE),nSpacing);
crlb_all = zeros(10,length(nTE),nSpacing);
% fim_all = zeros(10,10,length(TR_r),nSpacing);
% crlb_all = zeros(10,length(TR_r),nSpacing);
sigma = 0.5;
for kte = 1:length(nTE)
% for kte = 1:length(TR_r)
    for kSpacing = 1:nSpacing
        t = te_r(kte,kSpacing).te;
        [fim, crlb] = mwi_crlb_3cc_mwf(t,sigma,S0,mwf,icvf,r2s_mw,r2s_iw,r2s_ew,w_mb,w_ib,w_eb,phi);
        fim_all(:,:,kte,kSpacing) = fim;
        crlb_all(:,kte,kSpacing) = crlb;
    end
end
figure;imagesc(sqrt(squeeze(crlb_all(2,:,:))),[0 0.12]);
% yticklabels(ytickL);
% yticks(TR_r);