clear

syms R1IEW kIEWM M0EW M0IW M0MW TR FA R1M rhoMW

M0M     = M0MW ./ rhoMW;
M0IEW   = (M0EW + M0IW);

kMIEW = kIEWM * M0IEW/M0M;

Lambda_L = [-R1M-kMIEW  , kIEWM         ;
            kMIEW      	, -R1IEW-kIEWM  ];
        
Xi_L = expm(TR*Lambda_L);

% Use same notation as paper
C_L = [R1M*M0M; R1IEW*M0IEW];
I = eye(length(C_L));
    
% Use same notation as paper
T = cos(FA)*I;

% Spencer and Fishbein J. Magn. Reson:142(2000)120-135
steadystate_pool = sin(FA)* ((I-Xi_L*T) \ (Xi_L - I)) * (Lambda_L \ C_L);

SMW = steadystate_pool(1) * rhoMW;	        % MW = Myelin volume * MW density
SIW = steadystate_pool(2) * M0IW/(M0EW + M0IW);     % IW = IEW * v_ic
SEW = steadystate_pool(2) * M0EW/(M0EW + M0IW);     % EW = IEW * (1-v_ic)

% DS/DM0MW
dSMWdM0MW = simplify(diff(SMW,M0MW));
dSIWdM0MW = simplify(diff(SIW,M0MW));
dSEWdM0MW = simplify(diff(SEW,M0MW));
% DS/DM0IW
dSMWdM0IW = simplify(diff(SMW,M0IW));
dSIWdM0IW = simplify(diff(SIW,M0IW));
dSEWdM0IW = simplify(diff(SEW,M0IW));
% DS/DM0EW
dSMWdM0EW = simplify(diff(SMW,M0EW));
dSIWdM0EW = simplify(diff(SIW,M0EW));
dSEWdM0EW = simplify(diff(SEW,M0EW));
% DS/DkIEWM
dSMWdkIEWM = simplify(diff(SMW,kIEWM));
dSIWdkIEWM = simplify(diff(SIW,kIEWM));
dSEWdkIEWM = simplify(diff(SEW,kIEWM));
% DS/DR1IEW
dSMWdR1IEW = simplify(diff(SMW,R1IEW));
dSIWdR1IEW = simplify(diff(SIW,R1IEW));
dSEWdR1IEW = simplify(diff(SEW,R1IEW));
% DS/DR1IEW
dSMWdR1M = simplify(diff(SMW,R1M));
dSIWdR1M = simplify(diff(SIW,R1M));
dSEWdR1M = simplify(diff(SEW,R1M));

% 
% ss = ones(1,length(C_L)) * ss_pool;
% 
% dssdR1f = diff(ss,R1f);
%%
R1M_in = 4.27;
R1IEW_in = 1;
kIEWM_in = 6;
M0MW_in = 0.1;
M0IW_in = 0.6;
M0EW_in = 0.3;
rhoMW_in = 0.5;
TR_in = 50e-3;
FA_in = deg2rad(10);


tic
double(subs(dSMWdM0MW,[R1M R1IEW kIEWM M0MW M0IW M0EW rhoMW TR FA],[R1M_in R1IEW_in kIEWM_in M0MW_in M0IW_in M0EW_in rhoMW_in TR_in FA_in]))
toc
tic
mwi_ssSPGR_2T1_dSMWdM0MW(R1M_in,R1IEW_in,kIEWM_in,M0MW_in,M0IW_in,M0EW_in,rhoMW_in,TR_in,FA_in)
toc

tic
double(subs(dSIWdM0MW,[R1M R1IEW kIEWM M0MW M0IW M0EW rhoMW TR FA],[R1M_in R1IEW_in kIEWM_in M0MW_in M0IW_in M0EW_in rhoMW_in TR_in FA_in]))
toc
tic
mwi_ssSPGR_2T1_dSIWdM0MW(R1M_in,R1IEW_in,kIEWM_in,M0MW_in,M0IW_in,M0EW_in,rhoMW_in,TR_in,FA_in)
toc


% 
% dx = 1e-6;
% (double(subs(ss,[R1r R1f krf M0f M0r tr fa],[4.27 1 4 0.8 0.2 50e-3 5])) - ...
%     double(subs(ss,[R1r R1f krf M0f M0r tr fa],[4.27 1+dx 4 0.8 0.2 50e-3 5]))) ./ dx_in