%% [ss,ss_pool] = mwi_model_ssSPGR_3T1_X(fa,tr,...
%                                       M0a,M0b,M0c,...
%                                       t1a,t1b,t1c,...
%                                       kba,kbc,kca,kac)
%
% Input
% --------------
% fa            : flip angles (degree, sould be b1 corrected)
% tr            : repetition time (s or ms, must be same unit as T1s)
% M0a           : Compartment A signal
% M0b           : Compartment B signal
% M0c           : Compartment C signal
% t1a           : Compartment A T1 (s or ms, must be same unit as tr)
% t1b           : Compartment B T1 (s or ms, must be same unit as tr)
% t1c           : Compartment C T1 (s or ms, must be same unit as tr)
% kba           : exchange rate from B to A (s-1 or ms-1, must be same unite as tr)
% kbc           : exchange rate from B to C (s-1 or ms-1, must be same unite as tr)
% kca           : exchange rate from C to A (s-1 or ms-1, must be same unite as tr)
% kac           : exchange rate from A to C (s-1 or ms-1, must be same unite as tr)
%
% Output
% --------------
% ss            : total steady-state signal
% ss_pool       : steady-state signal coming from each pool, [A,B,C]
%
% Description: 3-pool exchange model for SPGR steady state
% ref: Dortch et al. 2013 MRM
% A system with 3 compartment can be in two forms:
%     A
%    / \
%   B - C (cyclic)    or A - B - C (linear)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 Feb 2019
% Date last modified: 
%
%% main
function [ss,ss_pool] = mwi_model_ssSPGR_3T1(fa,tr,...
                                      M0a,M0b,M0c,...
                                      t1a,t1b,t1c,...
                                      kba,kbc,kca,kac)
% set default values
if nargin < 9
    kba=0;
    kbc=0;
end
if nargin < 11
    kac=0;  
    kca=0; % A - B - C
end

% no. of flip angles
nfa = length(fa);

% derive the rest of the exchange rates
kab	= (kba + kbc) * M0b/M0a - kac;
kcb	= (kba + kbc) * M0b/M0c - kca;

% convert T1s to R1s for simplicity
R1a    = 1/t1a;
R1b    = 1/t1b;
R1c    = 1/t1c;

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [-R1a-kab-kac   , kba            , kca       	;
            kab            , -R1b-kba-kbc	, kcb          	;
            kac            , kbc            , -R1c-kcb-kca	];
        
Xi_L = expm(tr*Lambda_L);

% Use same notation as paper
C_L = [R1a*M0a; R1b*M0b; R1c*M0c];
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
