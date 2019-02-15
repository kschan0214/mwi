%%  [ss,ss_pool] = mwi_model_ssSPGR_4T1_linear(fa,tr,...
%                                       M0m,M0mw,M0iew,M0nm,...
%                                       t1m,t1mw,t1iew,t1nm,...
%                                       kmwm,kmwiew,kiewnm,...
%                                       b1)
%
% Input
% --------------
% fa            : flip angles (degree, sould be b1 corrected)
% tr            : repetition time
% M0m           : Myelin signal
% M0mw          : Myelin water signal
% M0iew         : Intra-axonal+extracellular water signal
% M0nm          : non-myelin macromolecules signal
% t1m           : myelin T1
% t1mw          : myelin water T1
% t1iew         : intra+extra-cellular water T1
% t1nm          : non-myelin macromolecules T1
% b1            : factor of B1+ inhomogeneity: true flip angle=1
%
% Output
% --------------
% ss            : total steady-state signal
% ss_pool       : steady-state signal coming from each pool, [m,mw,iew,nm]
%
% Description: Model based on Barta et al. J. Magn. Reson:259(2015)56-67
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 Feb 2019
% Date last modified: 
%
%% main
function [ss,ss_pool] = mwi_model_ssSPGR_4T1_linear(fa,tr,...
                                      M0m,M0mw,M0iew,M0nm,...
                                      t1m,t1mw,t1iew,t1nm,...
                                      kmwm,kmwiew,kiewnm)

% no. of flip angles
nfa = length(fa);

% derive the rest of the exchange rates
kmmw    = (kmwm     +   kmwiew) * M0mw /M0m;
kiewmw  = (kmwm     +   kmwiew) * M0mw /M0iew - kiewnm;
knmiew  = (kiewmw   +   kiewnm) * M0iew/M0nm;

% convert T1s to R1s for simplicity
R1m     = 1/t1m;
R1mw    = 1/t1mw;
R1iew   = 1/t1iew;
R1nm    = 1/t1nm;

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [-R1m-kmmw   , kmwm              , 0                     , 0     ;
            kmmw        , -R1mw-kmwm-kmwiew , kiewmw                , 0     ;
            0           , kmwiew            , -R1iew-kiewmw-kiewnm  , knmiew;
            0           , 0                 , kiewnm                , -R1nm-knmiew];
        
Xi_L = expm(tr*Lambda_L);

% Use same notation as paper
C_L = [R1m*M0m; R1mw*M0mw; R1iew*M0iew; R1nm*M0nm];
I = eye(length(C_L));

ss_pool = zeros(length(C_L),nfa);
ss = zeros(1,nfa);
for kfa = 1:length(fa)
    
    % Use same notation as paper
    T = cosd(fa(kfa))*I;
    
    % Spencer and Fishbein J. Magn. Reson:142(2000)120-135
    ss_pool(:,kfa) = sind(fa(kfa))* ((I-Xi_L*T) \ (Xi_L - I)) * (Lambda_L \ C_L);
    
    ss(kfa) = ones(1,length(C_L)) * ss_pool(:,kfa);

end

end
