function T3D = precomputeT_epg(phi,theta)
%   [F0,Fn,Zn,F] = EPG_GRE(theta,phi,TR,T1,T2,varargin)
%
%   Single pool EPG (classic version) for gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               phi:        phase per pulse. This function can hence be
%                           used to simulate RF spoiling or balanced
%                           sequences depending on how phase is cycled
%                           see function RF_phase_cycle()
%               TR:         repetition time, ms
%               T1:         T1, ms
%               T2:         T2, ms
%
%   optional arguments (use string then value as next argument)
%
%               kmax:       maximum EPG order to include. Can be used to
%                           accelerate calculation. 
%                           Setting kmax=inf ensures ALL pathways are
%                           computed
%              diff:        structure with fields:
%                           G    - Gradient amplitude(s)
%                           tau  - Gradient durations(s)
%                           D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%
%               prep:       can be used to simulate prep pulse prior to
%                           gradient echo. Assume all transverse
%                           magnetization after prep is spoiled.
%                           structure with fields:
%                           flip    -   flip angle, rad
%                           t_delay -   time delay, ms
%
%   Outputs:                
%               F0:         signal (F0 state) directly after each
%                           excitation
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0 F0* Z0 F1 F-1* Z1 F2 F-2* Z2 ...] etc
%
%
%   Shaihan Malik 2017-07-20

%%
%%% The maximum order varies through the sequence. This can be used to speed up the calculation    
np = length(phi);
% if not defined, assume want max
if ~exist('kmax','var')
    kmax = np - 1;
end

if isinf(kmax)
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = np - 1; % this is maximum value
else
    allpathways = false;
end

%%% Variable pathways
if allpathways
    kmax_per_pulse = 0:kmax;
else
    kmax_per_pulse = [1:ceil(np/2) (floor(np/2)):-1:1];
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
     
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end

%%% Number of states is 6x(kmax +1) -- +1 for the zero order
N=3*(kmax+1);

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

% store the indices of the top 3x3 corner, this helps build_T
i1 = [];
for ii=1:3
    i1 = cat(2,i1,sub2ind(size(T),1:3,ii*ones(1,3)));
end

%% Main body of gradient echo sequence, loop over TRs 

for jj=1:np 
    %%% RF transition matrix
    A = RF_rot(theta,phi(jj));
    
    %%% Replicate A to make large transition matrix
    build_T(A);
    
    T3D{jj}=T;
    
end

    %%% NORMAL EPG transition matrix as per Weigel et al JMR 2010 276-285 
    function Tap = RF_rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
    end

    function build_T(AA)
        ksft = 3*(3*(kmax+1)+1);
        for i2=1:9
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end
