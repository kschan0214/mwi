classdef HCFM
% This class implements the equations in 
% Wharton S, Bowtell R. NeuroImage. 2013;83:1011-1023. doi:10.1016/j.neuroimage.2013.07.054
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 19 Sep 2023
% Date modified:
%
    properties (Constant)
        gyro = 2*pi*42.57747892;
    end

    properties (GetAccess = public, SetAccess = protected)
        t;
        B0 = 3; %T
        normFactor;

    end

    methods

        function obj    = HCFM(t,B0)
            if nargin == 1 || ~isempty(t)
                obj.t   = t(:);
            end
            if nargin == 2 || ~isempty(B0)
                obj.B0          = B0;
                obj.normFactor  = obj.gyro*obj.B0;
            end

        end

        function t      = TransitionTimeExtraaxonal(obj,x_i,x_a,g,theta)
        %
        % Input
        % --------------
        % x_i           : isotropic susceptibility, in ppm
        % x_a           : anisotropic susceptibility, in ppm
        % g             : g-ratio, [0,1]
        % theta         : angle between B0 and fibre direction, radian
        % b0            : field strength, T
        %
        % Output
        % --------------
        % t             : transition time of dephasing from quadratic regime to
        %                 linear regime
        %
        % Description:
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 12 July 2019
        % Date modified:
        %
            % Eq. A9
            x_d = obj.EffectiveSusceptibility(x_i,x_a,g);
            
            % Eq. A8
            % t = 3 ./ (abs(x_d)*obj.gyro*obj.B0.*sin(theta).^2); % this is the solution provided in the paper
            % the true transition time should be:
            t = 4 ./ (abs(x_d)*obj.gyro*obj.B0.*sin(theta).^2); 
        
        end

        function D_E    = DephasingExtraaxonal(obj,fvf,g,x_i,x_a,theta)
        %
        % Input
        % --------------
        % fvf           : fibre volume fraction, [0 1]
        % g             : g-ratio, [0 1]
        % x_i           : isotropic susceptibility, ppm
        % x_a           : anisotropic susceptibility, ppm
        % theta         : angle between B0 and fibre direction, radian
        %
        % Output
        % --------------
        % D_E           : decay weight at each time point
        %
        % Description:
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 12 July 2019
        % Date modified:
        %
            % Eq.A8 Wharton.NI2013
            alpha = obj.TransitionTimeExtraaxonal(x_i,x_a,g,theta);
            % Eq.A9 Wharton.NI2013
            x_d   = obj.EffectiveSusceptibility(x_i,x_a,g);
            
            D_E = zeros(size(obj.t));
            for kt=1:length(obj.t)
                if obj.t(kt) < alpha
                    D_E(kt) = (fvf/16)*(abs(x_d)*obj.gyro*obj.B0*sin(theta).^2*obj.t(kt))^2;
                else
                    D_E(kt) = (fvf/2)*(abs(x_d)*obj.gyro*obj.B0*sin(theta).^2)*...
                        (obj.t(kt) - 2./(abs(x_d)*obj.gyro*obj.B0*sin(theta).^2));
                end
            end
        end

        function freq   = FrequencyMyelin(obj,x_i,x_a,g,theta,E)
        %
        % Input
        % --------------
        % x_i           : isotropic susceptibility, ppm
        % x_a           : anisotropic susceptibility, ppm
        % g             : g-ratio, [0,1]
        % theta         : angle between B0 and fibre direction, radian
        % E             : exchange, ppm
        %
        % Output
        % --------------
        % freq_myelin   : frequency shift in myelin, radian
        %
        % Description:
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 12 July 2019
        % Date modified:
        %
            % Eq. [A15]
            c1 = obj.C1(g);
            
            % Eq.[5]
            freq = ((x_i./2).*(2/3 - sin(theta).^2) + (x_a./2).*(c1.*sin(theta).^2 - 1/3) + E)*obj.gyro*obj.B0;

        end

        function freq   = FrequencyAxon(obj,x_a,g,theta)
        %
        % Input
        % --------------
        % x_a           : anisotropic susceptibility, ppm
        % g             : g-ratio, [0,1]
        % theta         : angle between B0 and fibre direction, radian
        % b0            : field strength, in T
        %
        % Output
        % --------------
        % freq          : frequency shift in axon, in radian
        %
        % Description:
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 12 July 2019
        % Date modified:
        %
        %
            % Eq.[6]
            freq = (3*x_a./4) .* sin(theta).^2 .* log(1./g) * obj.gyro*obj.B0;

        end

    end

    methods (Static)

        function x_d    = EffectiveSusceptibility(x_i,x_a,g)
        %
        % Input
        % --------------
        % x_i           : isotropic susceptibility, ppm
        % x_a           : anisotropic susceptibility, ppm
        % g             : g-ratio, [0 1]
        % 
        % Output
        % --------------
        % x_d           : effective susceptibility, in ppm
        %
        % Description: Eq.A9, Wharton and Bowtell 2013 NI
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 12 July 2019
        % Date modified: 13 Sep 2023
        %
            x_d = (x_i + x_a/4).*(1-g.^2);

        end

        function c1     = C1(g)
        %
        % Input
        % --------------
        % g             : g-ratio
        %
        % Output
        % --------------
        % c1            : c1 of equation A15
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 12 July 2019
        % Date last modified:
        %
        %
            c1          = 1/4 - (3/2)*((g.^2)./(1-g.^2)).*log(1./(g.^2));
            c1(g==1)    = 1/4 - 3/2;
        end

        function [v_myelin,v_axon,v_ec] = VolumeFraction(fvf,g)
        %
        % Input
        % --------------
        % fvf           : fibre volume fraction, [0,1]
        % g             : g-ratio, [0,1]
        %
        % Output
        % --------------
        % v_myelin      : volume fraction of myelin
        % v_axon        : volume fraction of intra-axonal space
        % v_ec          : volume fraction of extra-axonal space
        %
        % Description:
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 12 July 2019
        % Date modified:
        %
            v_myelin	= fvf .* (1-g.^2);
            v_axon  	= fvf .* g.^2;
            v_ec        = 1-fvf;
        
        end
    end

end