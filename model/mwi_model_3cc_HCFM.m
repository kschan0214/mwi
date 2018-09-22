%% s = mwi_model_3cc_HCFM(te, fvf, g, rho_m, e, ...
%%                        t2_m, t2_n, chi_I, chi_A, theta,...
%%                        freq_bkg, pini, s0, b0)
%
% Input
% --------------
% te            : echo times, in second
% fvf           : fibre volume fraction, [0,1]
% g             : g-ratio = radius_inner/radius_outer, [0,1]
% rho_m         : relative concentration of myelin water to the
%                 extra-cellular and intra-axonal water
% e             : orientation-independent parameter decribing the
%                 contribution of chemical exchange within the myelin sheath 
%                 to the frequency offset in this region, in ppm
% t2_m          : T2 of myelin water, in second
% t2_n          : T2 of non-myelin water, in second
% chi_I         : isotropic magnetic susceptibility of myelin sheath, in
%                 ppm
% chi_A         : anisotropic magnetic susceptibility along principal axis,
%                 of anisotropic susceptibility tensor, in ppm
% theta         : fibre orientation w.r.t main magnetic field
% s0            : average signal strength
% f_bkg         : background field inhomogeneities, in Hz
% pini          : inital signal phase introduced by B1, [-pi,pi)
% b0            : magnetic field strength, in Tesla
%
% Output
% --------------
% s             : complex-valued signal with the same size as te
%
% Description: return complex-valued MR signals with contribution of myelin
% water, intra-axonal water and extra-cellular water. For more information,
% please check Wharton and Bowtell, NeuroImage:83 1011-1023(2013)
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 23 May 2018
% Date last modified:
%
%
function s = mwi_model_3cc_HCFM(te, fvf, g, rho_m, e, ...
                                t2_m, t2_n, chi_I, chi_A, theta,...
                                freq_bkg, pini, s0, b0)

if nargin < 14
    b0 = 3;             % T
end

if nargin < 13
    s0 = 1;             % a.u.
end

if nargin < 12
    pini = 0;           % rad
end

if nargin < 11
    freq_bkg = 0;       % Hz
end

gyro = 42.57747892 * 2*pi;     % 10^6 rad s^-1 T^-1


te = te(:);             % s
nt = length(te);

% effective susceptibility in extra-cellular water, Eq.[A9]
chi_D = (chi_I + chi_A/4)*(1-g^2);

% time of dephasing transition from quadratic to linear regimes, Eq.[A8]
alpha = 3/(abs(chi_D)*gyro*b0*sind(theta)^2);

% coefficients, Eq.[A15]
c1 = 1/4 - (3/2)*((g^2)/(1-g^2))*log(1/g);

% myelin water average frequency shift, Eq.[5]
w_m = ((chi_I/2) * (2/3-sind(theta)^2) + ...
       (chi_A/2) * (c1*sind(theta)^2-1/3) + e) * gyro*b0;
   
% intra-axonal water average frequency shift, Eq[6]
w_a = (3*chi_A*sind(theta).^2)/4 * log(1/g)    * gyro*b0;

% assuming extracellular water is on resonance
w_e = 0;            

% angilar frequency of background field
w_bkg = freq_bkg * 2*pi;
                                   
% linear variables
v_m = fvf * (1-g^2);    % Eq.[A1]
v_a = fvf * (g^2);      % Eq.[A2]
v_e = (1-fvf);          % Eq.[A3]

% ommit decay effect of signal dephasing of myelin water
d_m = 0;  

s_m = zeros(nt,1, 'like', te);
s_a = zeros(nt,1, 'like', te);
s_e = zeros(nt,1, 'like', te);
for kt = 1:nt
    
    % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
    constant = abs(chi_D)*gyro*b0*sind(theta)^2;
    if te(kt) < alpha
        d_e = (fvf/16)*(constant*te(kt))^2;
    else
        d_e = (fvf/2) *constant*(te(kt)-2/constant);
    end
    
    % myelin water signal, Eq.[2]
    s_m(kt) = v_m*rho_m * exp(-te(kt)/t2_m) * exp(-d_m) ...  % decay
                        * exp(1i*(w_m+w_bkg)*te(kt));   % phase change
    
    % intra-axonal water, Eq.[3]
    s_a(kt) = v_a       * exp(-te(kt)/t2_n) ...              % decay
                        * exp(1i*(w_a+w_bkg)*te(kt));   % phase change
                   
    % extra-cellular water, Eq.[4]
    s_e(kt) = v_e       * exp(-te(kt)/t2_n) * exp(-d_e) ...  % decay
                        * exp(1i*(w_e+w_bkg)*te(kt));   % phase change
                   
end

s = s0 * (s_m + s_a + s_e) * exp(-1i*pini);  % Eq[1] + initial phase due to B1



% frequency dephasing
% psi = mwi_psi_t2sw([0, 0, 0], [f_m+f_bkg, f_a+f_bkg, f_e+f_bkg], te);
% 
% s = psi * [v_m*pd_m; v_a; v_e] * exp(-1i*pini);



end