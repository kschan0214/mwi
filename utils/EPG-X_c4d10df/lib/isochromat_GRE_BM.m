function [s,mxy] = isochromat_GRE_BM(theta,phi,TR,T1x,T2x,f,ka,Niso,varargin)
%   [F0,Fn,Zn,F] = isochromat_GRE_BM(theta,phi,TR,T1x,T2x,f,ka,Niso,varargin)
%
%   Isochromat based simualtion of Bloch McConnell coupled systems w/ gradient echo sequences
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
%               Niso:       Number of isochromats for simulation
%
%
%   Optional:                
%               delta:      frequency offset of pool b, kHz 
%
%   Shaihan Malik 2017-09-04
%                 2017-10-06: rewrite in M+,M- basis  

%% Extra variables

for ii=1:length(varargin)
    
    
    % frequency offset of pool b, for phase gain during evolution period
    if strcmpi(varargin{ii},'delta')
        delta = 2*pi*varargin{ii+1};
    end
    
    % phase range from unbalanced gradient, default is 2pi
    if strcmpi(varargin{ii},'phase_range')
        phase_range = varargin{ii+1};
    end
    
end

if ~exist('delta','var')
    delta=0;
end

%% Set up variables


%%% Fix number of TR periods
Npulse = length(theta);

%%% Isochromat phase distribution
if ~exist('phase_range','var')
    phase_range = 2*pi;
end
psi = phase_range*(0:fix(Niso)-1)/fix(Niso);

%%% Number of variables (each isochromat has Mxa,Mya,Mza,Mxb,Myb,Mzb)
N = 6*Niso;

%% Dependent variables for exchange case
M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;

R1a = 1/T1x(1);
R1b = 1/T1x(2);
R2a = 1/T2x(1);
R2b = 1/T2x(2);

%%% handle no exchange case
if (f==0)||(ka==0)
    ka = 0;
    kb = 0;
    M0b = 0;
    R1b = 1e-3;%<- doesn't matter what it is, just avoid singular matrices
end


%% Set up matrices for Relaxation and Exchange

%%% Relaxation-exchange matrix
La = diag(-[R2a R2a R1a]);
Lb = diag(-[R2b R2b R1b]);
Ka = diag([ka ka ka]);
Kb = diag([kb kb kb]);
A = [[La-Ka Kb];[Ka Lb-Kb]];
C = [0 0 R1a*M0a 0 0 R1b*M0b]';

Xi = expm(TR*A); %<-- operator for time period TR


% For inhomogeneous solution also need:
Zoff = (Xi - eye(6))*(A\C);
Zoff = repmat(Zoff(:),[Niso 1]); % replicate for all isochromats

% 6/10/17: make Xi include dephasing and delta_b
Xi={};
w = psi/(TR);
for jj=1:Niso
    W = diag(1i*[-w(jj) w(jj) 0 -(w(jj)+delta) w(jj)+delta 0]);
    AA = A + W;
    Xi{jj} = expm(AA*TR);
end
Xi = blkdiag(Xi{:});
Xi = sparse(Xi);  

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

%%% Matrix to hold all results
M = zeros([N Npulse]);

% Initialize Mz=1 for all isochromats
M0 = zeros([N 1]);
M0(3:6:N)=M0a;
M0(6:6:N)=M0b;
M(:,1) = M0;

%% Main body of gradient echo sequence, loop over TRs 

% Loop over pulses: flip, record FID then dephase
for jj=1:Npulse
    
    % Get matrix ready for flip
    build_T_matrix_sub_implicit(RF_rot(theta(jj),phi(jj)));
    
    % Flip
    M(:,jj) = T*M(:,jj);
    
    % Dephase and T1 recovery
    if jj<Npulse
      
        % Evolution due to gradients, relaxation and exchange all together
        M(:,jj+1) = Xi*M(:,jj)+Zoff;
    end
     
end


% Now generate signal and demodulate it
% mxya = M(1:6:end,:) + 1i*M(2:6:end,:);
% mxyb = M(4:6:end,:) + 1i*M(5:6:end,:);
mxya = M(1:6:end,:);
mxyb = M(4:6:end,:);

% Get signal from mean
sa = 1i*mean(mxya,1);
sb = 1i*mean(mxyb,1);
% demodulate this
sa = sa .* exp(-1i*phi(1:Npulse));
sb = sb .* exp(-1i*phi(1:Npulse));

s=sa+sb;%return total signal
mxy=cat(3,mxya,mxyb);


    function build_T_matrix_sub_implicit(AA)
        %%% This function operates on the existing T matrix, rather than
        %%% re-declare it each time. This is much faster 
        ix = 1:(N/3);
        T(3*N*(ix-1)+3*(ix-1)+1)=AA(1);
        T(3*N*(ix-1)+3*(ix-1)+2)=AA(2);
        T(3*N*(ix-1)+3*(ix-1)+3)=AA(3);
        T(3*N*ix-2*N+3*(ix-1)+1)=AA(4);
        T(3*N*ix-2*N+3*(ix-1)+2)=AA(5);
        T(3*N*ix-2*N+3*(ix-1)+3)=AA(6);
        T(3*N*ix-N+3*(ix-1)+1)=AA(7);
        T(3*N*ix-N+3*(ix-1)+2)=AA(8);
        T(3*N*ix-N+3*(ix-1)+3)=AA(9);
    end

    % Rotation matrix function
    function R = rotmat(u)
        
        % check zero input
        if any(u)
            th=norm(u);
            u=u/th;
            ct=cos(th);
            st=sin(th);
            R = [[ct + u(1)^2*(1-ct) u(1)*u(2)*(1-ct)-u(3)*st u(1)*u(3)*(1-ct)+u(2)*st];...
                [u(2)*u(1)*(1-ct)+u(3)*st ct+u(2)^2*(1-ct) u(2)*u(3)*(1-ct)-u(1)*st];
                [u(3)*u(1)*(1-ct)-u(2)*st u(2)*u(3)*(1-ct)+u(1)*st] ct+u(3)^2*(1-ct)];
            
        else
            R=eye(3);
        end
        
    end

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
        %Tap = kron(eye(2),Tap);
    end
end
