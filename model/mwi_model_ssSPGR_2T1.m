%% [ss,ss_pool] = mwi_model_ssSPGR_2T1(fa,tr,...
%                                       M0r,M0f,...
%                                       t1r,t1f,...
%                                       krf,kfr)
%
% Input
% --------------
% fa            : flip angles (degree, sould be b1 corrected)
% tr            : repetition time
% M0mw          : Myelin water signal
% M0iw          : Intra-axonal water signal
% M0ew          : extracellular water signal
% t1mw          : myelin water T1
% t1iw          : intra-axonal water T1
% t1ew          : extracellular water T1
% b1            : factor of B1+ inhomogeneity: true flip angle=1
%
% Output
% --------------
% ss            : total steady-state signal
% ss_pool       : steady-state signal coming from each pool, 
%                 1st dimension - [r,f]; 2nd dim - flip angle
%
% Description: Model based on Ou and Gochberg 2008 MRM
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 Feb 2019
% Date modified: 17 Nov 2020, bug fix when M0r or M0f = 0
%
%% main
function [ss,ss_pool] = mwi_model_ssSPGR_2T1(fa,tr,M0r,M0f,t1r,t1f,krf,kfr)

% no. of flip angles
nfa = length(fa);

% derive the rest of the exchange rates
if isempty(krf)
    if M0r ~= 0
        krf = kfr * M0f/M0r;
    else
        krf = 0;
    end
        
else
    if M0f ~= 0
        kfr   = krf * M0r/M0f;
    else
        kfr = 0;
    end
end

% convert T1s to R1s for simplicity
R1r    = 1/t1r;
R1f    = 1/t1f;

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [-R1r-krf    , kfr       ;
            krf      	, -R1f-kfr  ];
        
Xi_L = expm(tr*Lambda_L);

% Use same notation as paper
C_L = [R1r*M0r; R1f*M0f];
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
