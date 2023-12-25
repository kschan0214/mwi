%% Initial guesses for MWI @ 3T, In vivo

%%%%%%%%%% Signal intensity %%%%%%%%%%
% Myelin water fraction, if not provided then set to 10%
if ~exist('mwf0','var')
    mwf0 = 0.1; 
end  
% Intra-axonal volume fraction, if not provided then set to 70%
if ~exist('v_ic','var')
    v_ic = 0.7;
end
% derived parameter
iwf0 = (1-mwf0)*v_ic;	% Intra-axonal water fraction
ewf0 = 1-mwf0-iwf0;  	% Extra-axonal water fraction

%%%%%%%%%% R2* %%%%%%%%%%
% MW R2*, corresponding to T2* ranging from [3.3, 20]ms, starting at 10ms 
r2smwlb = 50;	% lower bound
r2smwub = 300;  % upper bound
r2smw0  = 100;  % Starting point

% IEW R2*
if ~exist('r2siew0','var') || isnan(r2siew0) || isinf(r2siew0)
    r2siew0 = 18.5; % 3T
end
% IW R2*, corresponding to T2* ranging from [20, 167]ms, starting at 63ms by default or previously estimated value 
r2siwlb = 6;	
r2siwub = 50;   
r2siw0  = min(r2siew0-2.5,r2siwub); r2siw0 = max(r2siw0,r2siwlb); % make sure initial point is within range
% EW R2*, corresponding to T2* ranging from [20, 167]ms, starting at 48ms by default or previously estimated value 
r2sewlb = 6;	
r2sewub = 50;   
r2sew0  = min(r2siew0+2.5,r2sewub); r2sew0 = max(r2sew0,r2sewlb); % make sure initial point is within range

%%%%%%%%%% T1 %%%%%%%%%%
% MW T1
t1slb   = 50e-3;  	
t1sub   = 500e-3;           
t1s0    = 234e-3;

% IEW T1
if ~exist('t1iew0','var') || isnan(t1iew0) || isinf(t1iew0)
    t1iew0 = 1; % 3T
end
t1llb   = 500e-3; 	
t1lub   = min(t1iew0+1,5);	
t1l0    = min(t1iew0,t1lub); t1l0   = max(t1l0,t1llb); % make sure initial point is within range

% exchange rate from IEW T1 to MW T1, s^-1
klslb   = 0;      	
klsub   = 20;               
kls0    = 0;       	 
