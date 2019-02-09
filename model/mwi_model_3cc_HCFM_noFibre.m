function s = mwi_model_3cc_HCFM_noFibre(te,fvf,g,t2s_mw,t2_ow,w_mw,s0,w_bkg,x_i,x_a,pini,rho_mw,E,b0)

if nargin < 14;	b0 = 3;         end
if nargin < 13; E = 0.02;       end
if nargin < 12; rho_mw = 0.7; 	end
if nargin < 11; pini = 0;       end
if nargin < 10; x_a = 0.1;      end
if nargin < 9;  x_i = 0.1;      end
if nargin < 8;  w_bkg = 0;      end
if nargin < 7;  s0 = 1;         end

nt = length(te);

gyro = 42.57747892 * 2*pi;     % 10^6 rad s^-1 T^-1

% coefficients, Eq.[A15]
c1 = 1/4 - (3/2)*((g^2)/(1-g^2))*log(1/g);

% effective susceptibility in extra-cellular water, Eq.[A9]
x_D = (x_i + x_a/4)*(1-g^2);

% sin(theta)^2 derived using myelin frequency
sine_sq_theta = (6*w_mw - 2*x_i + (3*E-1)*x_a) / (3*gyro*b0*(c1*x_a-x_i));

w_iw = log(1/g)*gyro*b0*(3*x_a*sine_sq_theta)/4;

% time of dephasing transition from quadratic to linear regimes, Eq.[A8]
alpha = 3/(abs(x_D)*gyro*b0*sine_sq_theta);

% linear variables
v_m = fvf * (1-g^2);    % Eq.[A1]
v_a = fvf * (g^2);      % Eq.[A2]
v_e = (1-fvf);          % Eq.[A3]

s_mw = zeros(nt,1, 'like', te);
s_iw = zeros(nt,1, 'like', te);
s_ew = zeros(nt,1, 'like', te);
for kt = 1:nt
    
    % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
    constant = abs(x_D)*gyro*b0*sine_sq_theta;
    if te(kt) < alpha
        d_e = (fvf/16)*(constant*te(kt))^2;
    else
        d_e = (fvf/2) *constant*(te(kt)-2/constant);
    end
    
    % myelin water signal, Eq.[2]
    s_mw(kt) = v_m*rho_mw * exp(-te(kt)/t2s_mw) ...  % MW T2* decay
                          * exp(1i*(w_mw+w_bkg)*te(kt));   % phase change
    
    % intra-axonal water, Eq.[3]
    s_iw(kt) = v_a        * exp(-te(kt)/t2_ow) ...              % T2 decay
                          * exp(1i*(w_iw+w_bkg)*te(kt));   % phase change
                   
    % extra-cellular water, Eq.[4]
    s_ew(kt) = v_e        * exp(-te(kt)/t2_ow) * exp(-d_e) ...  % T2 decay
                          * exp(1i*(w_bkg)*te(kt));   % phase change
                   
end

s = s0 * (s_mw + s_iw + s_ew) * exp(-1i*pini);  % Eq[1] + initial phase due to B1

end