function  T3D = PrecomputeT(phi,B1)
%   [F0,Fn,Zn,F] = EPGX_GRE_BM(theta,phi,TR,T1x,T2x,f,ka,varargin)
%
%   EPG-X for Bloch McConnell coupled systems w/ gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               phi:        phase per pulse. This function can hence be
%                           used to simulate RF spoiling or balanced
%                           sequences depending on how phase is cycled
%                           see function RF_phase_cycle()
%               TR:         repetition time, ms
%               T1x:        [T1a T1b], ms
%               T2x:        [T2a T2b], ms
%               f:          fraction of compartment b
%               ka:         forward exchange rate from a->b (units ms^-1)
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
%           * Diffusion is same in both compartments, this is experimental *
%
%               prep:       can be used to simulate prep pulse prior to
%                           gradient echo. Assume all transverse
%                           magnetization after prep is spoiled.
%                           structure with fields:
%                           flip    -   flip angle, rad
%                           t_delay -   time delay, ms
%
%
%               delta:      frequency offset of pool b, kHz
%
%   Outputs:
%               F0:         signal (F0 state) directly after each
%                           excitation (saturation factor(s))
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0a F0a* Z0a F0b F0b* Z0b F1a F-1a* Z1a F1b F-1b* Z1b ...] etc
%
%
%   Shaihan Malik 2017-07-20


%% Extra variables

%%% The maximum order varies through the sequence. This can be used to speed up the calculation
np = length(phi);
% if not defined, assume want max
if ~exist('kmax','var')
    kmax = np - 1;
end


%%% Number of states is 6x(kmax +1) -- +1 for the zero order
N=6*(kmax+1);


%% Set up matrices for Relaxation and Exchange
    %%% Pre-allocate RF matrix
% T3D=zeros(N,N,np);
% T3D = sparse(T3D);
     T = zeros(N,N);
    T = sparse(T);
% T3D =repmat(T,[1 1 np]);
for hh=1:length(B1)
    %%% Pre-allocate RF matrix
    T = zeros(N,N);
    T = sparse(T);
    
    % store the indices of the top 6x6 corner, this helps build_T
    i1 = [];
    for ii=1:6
        i1 = cat(2,i1,sub2ind(size(T),1:6,ii*ones(1,6)));
    end
    
    
    %% Main body of gradient echo sequence, loop over TRs
    
    for jj=1:np
        %%% RF transition matrix
        A = RF_rot(B1(hh),phi(jj));
        
        %%% Replicate A to make large transition matrix
        build_T(A);
%         T3D(:,:,jj)=T;
        T3D{jj}=T;
        %%% Apply flip and store this: splitting these large matrix
        %%% multiplications into smaller ones might help
    end
    
%     save(['TmatrixLibrary/',num2str(round(B1(hh)))], 'T3D');
end
% save(['TmatrixLibrary/flipangleTable'],'B1');
%%% NORMAL EPG transition matrix but replicated twice
% As per Weigel et al JMR 2010 276-285
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
        Tap = kron(eye(2),Tap);
    end

% New version, rotates pool b differently
    function Tap = RF_rot_v2(a,p)
        Tap = zeros([6 6]);
        
        % pool a
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(7) = conj(Tap(2));
        Tap(8) = Tap(1);
        Tap(9) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(13) = -1i*exp(1i*p)*sin(a);
        Tap(14) = 1i*exp(-1i*p)*sin(a);
        Tap(15) = cos(a);
        
        % pool b
        a = a*offres;
        Tap(22) = cos(a/2).^2;
        Tap(23) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(24) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(28) = conj(Tap(23));
        Tap(29) = Tap(22);
        Tap(30) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(34) = -1i*exp(1i*p)*sin(a);
        Tap(35) = 1i*exp(-1i*p)*sin(a);
        Tap(36) = cos(a);
    end

    function build_T(AA)
        ksft = 6*(6*(kmax+1)+1);
        for i2=1:36
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end
