classdef GREMWI
% This class implements the equations in 
% Wharton S, Bowtell R. NeuroImage. 2013;83:1011-1023. doi:10.1016/j.neuroimage.2013.07.054
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 19 Sep 2023
% Date modified:
%
    properties (Constant)
        gyro = 42.57747892;
    end

    properties (GetAccess = public, SetAccess = protected)
        t;
        B0  = 3; % T
        normFactor;
        x_i; % ppm
        x_a; % ppm
        E  ; % ppm
        rho_mw;
        
    end

    methods

        function obj    = GREMWI(t,B0,x_i,x_a,E,rho_mw)
            if nargin == 1 || ~isempty(t)
                obj.t   = t(:);
            end
            if nargin == 2 || ~isempty(B0)
                obj.B0          = B0;
                obj.normFactor  = obj.gyro*obj.B0;
            end
            if nargin == 3 || ~isempty(x_i);    else; obj.x_i       = -0.1; end % ppm
            if nargin == 4 || ~isempty(x_a);    else; obj.x_a       = -0.1; end % ppm
            if nargin == 5 || ~isempty(E);      else; obj.E         = 0.02; end % ppm
            if nargin == 6 || ~isempty(rho_mw); else; obj.rho_mw    = 0.02; end % ppm

        end

        % Fitting functions
        function [E, J] = residuals(obj, pars, DIMWI)

            Amw = pars(1);    
            if DIMWI.isVic
                Aie = pars(2); 
                Aiw = Aie * DIMWI.icvf; Aew = Aie*(1-DIMWI.icvf);
                counter = 3;
            else
                Aiw = pars(2); Aew = pars(3); 
                counter = 4;
            end
            
            % T2*s
            t2smw = 1/x(counter);         counter = counter + 1;
            t2siw = 1/x(counter);         counter = counter + 1;
            if ~DIMWI.isR2sEW
                t2sew = 1/x(counter);   counter = counter + 1;
            else
                t2sew = t2siw;
            end
            
            % frequency shifts
            if ~DIMWI.isFreqMW
                freq_mwbg = x(counter)/(2*pi); counter = counter + 1;
            else
                freq_mwbg = 0;
            end
            if ~DIMWI.isFreqIW
                freq_iwbg = x(counter)/(2*pi); counter = counter + 1;
            else
                freq_iwbg = 0;
            end
            
            % external scanner effects
            if fitAlgor.numMagn==numel(te) % magnitude fitting
                fbg=0;        pini=0;
            else    % other fittings
                fbg=x(counter)/(2*pi);     pini=x(counter+1);
            end
            
            freq_mw = freq_mwbg - fbg;
            freq_iw = freq_iwbg - fbg;
            freq_ew = 0;

            shat    = zeros([numel(te), numel(DIMWI.theta)]);
            dshat   = zeros([numel(te), numel(DIMWI.theta)]);
            for ktheta = 1:length(DIMWI.theta)

                [shat(:,ktheta), dshat(:,ktheta)] = obj.Signal_T2sMWI_DIMWI(Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,theta,DIMWI); 

            end
            shat = sum(shat .* DIMWI.ff(:).');

            E = shat(:) - y(:);

            if nargout > 1
                J = dshat;
            end
        end
        
        % Signal generator
        function [s,ds] = Signal_T2sMWI(obj,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini)
        %
        % Input
        % --------------
        % Amw           : Myelin water signal
        % Aiw           : Axonal water signal
        % Aew           : extracellular water signal
        % t2smw         : myelin water T2*, in second
        % t2siw         : axonal water T2*, in second
        % t2sew         : extracellular water T2*, in second
        % freq_mw    	: myelin water frequency, in Hz
        % freq_iw   	: axonal water frequency, in Hz
        % freq_ew       : extracellular frequency, in Hz
        % fbg           : background field, in Hz
        % pini          : initial phase introduced by B1+ phase offset, in rad
        %
        % Output
        % --------------
        % S             : MGRE MWI signal
        %
        % Description: Complex-valued model fitted to complex-valued data used in
        % Nam et al. NeuroImages 2015 116:214-221
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 4 January 2018
        % Date modified: 16 August 2018
        % Date modified: 29 October 2019
        %

            te = obj.t;

            % combine compartmental frequency with background
            freq_mwbg = freq_mw + fbg;
            freq_iwbg = freq_iw + fbg;
            freq_ewbg = freq_ew + fbg;

            s_mw = obj.GRE_signal(Amw, -1/t2smw, freq_mwbg, te);
            s_iw = obj.GRE_signal(Aiw, -1/t2siw, freq_iwbg, te);
            s_ew = obj.GRE_signal(Aew, -1/t2sew, freq_ewbg, te);
    
            s = (s_mw + s_iw + s_ew) .* exp(1i*pini);
           
            % Jacobian
            if nargout == 2
                dsdAmw       = exp(te .* (-1/t2smw + 1i*2*pi*freq_mwbg)).* exp(1i*pini);
                dsdAiw       = exp(te .* (-1/t2siw + 1i*2*pi*freq_iwbg)).* exp(1i*pini);
                dsdAew       = exp(te .* (-1/t2sew + 1i*2*pi*freq_ewbg)).*exp(-d_e).* exp(1i*pini);
                dsdR2smw     = -te .* Amw * exp(te .* (-1/t2smw + 1i*2*pi*freq_mwbg)) .* exp(1i*pini);
                dsdR2siw     = -te .* Aiw * exp(te .* (-1/t2siw + 1i*2*pi*freq_iwbg)) .* exp(1i*pini);
                dsdR2sew     = -te .* Aew * exp(te .* (-1/t2sew + 1i*2*pi*freq_ewbg)).*exp(-d_e) .* exp(1i*pini);
                dsdfreq_mwbg = 1i*2*pi*te .* Amw * exp(te .* (-1/t2smw + 1i*2*pi*freq_mwbg)) .* exp(1i*pini);
                dsdfreq_iwbg = 1i*2*pi*te .* Aiw * exp(te .* (-1/t2siw + 1i*2*pi*freq_iwbg)) .* exp(1i*pini);
                dsdfreq_ewbg = 1i*2*pi*te .* Aew * exp(te .* (-1/t2sew + 1i*2*pi*freq_ewbg)).*exp(-d_e) .* exp(1i*pini);
                dsdpini      = 1i* ( Amw * exp(te .* (-1/t2smw + 1i*2*pi*freq_mwbg)) + ...
                                     Aiw * exp(te .* (-1/t2siw + 1i*2*pi*freq_iwbg)) + ...
                                     Aew * exp(te .* (-1/t2sew + 1i*2*pi*freq_ewbg)).*exp(-d_e) ).* exp(1i*pini);;

                ds = cat(1,dsdAmw(:),dsdAiw(:),dsdAew(:),dsdR2smw(:),dsdR2siw(:),dsdR2sew(:),...
                            dsdfreq_mwbg(:),dsdfreq_iwbg(:),dsdfreq_ewbg(:),dsdpini(:));
            end
        
        end

        function [s,ds] = Signal_T2sMWI_DIMWI(obj,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,theta,DIMWI)
        %
        % Input
        % --------------
        % Amw           : Myelin water signal
        % Aiw           : Axonal water signal
        % Aew           : extracellular water signal
        % t2smw         : myelin water T2*, in second
        % t2siw         : axonal water T2*, in second
        % t2sew         : extracellular water T2*, in second
        % freq_mw    	: myelin water frequency, in Hz
        % freq_iw   	: axonal water frequency, in Hz
        % freq_ew       : extracellular frequency, in Hz
        % fbg           : background field, in Hz
        % pini          : initial phase introduced by B1+ phase offset, in rad
        %
        % Output
        % --------------
        % S             : MGRE MWI signal
        %
        % Description: Complex-valued model fitted to complex-valued data used in
        % Nam et al. NeuroImages 2015 116:214-221
        %
        % Kwok-Shing Chan @ DCCN
        % kwokshing.chan@donders.ru.nl
        % Date created: 4 January 2018
        % Date modified: 16 August 2018
        % Date modified: 29 October 2019
        %

            te = obj.t;

            %%%%%%%%%% frequency shifts estimated using HCFM %%%%%%%%%%
            HCFM_obj = HCFM(obj.t, obj.B0);
            % derive g-ratio 
            g       = obj.gratio(abs(Amw)/obj.rho_mw, abs(Aiw));
            % derive fvf from signal intensity
            v_ic    = abs(Aiw) ./ (abs(Aiw) + abs(Aew));
            fvf     = v_ic ./ (g^2 - v_ic*g^2 + v_ic); 

            % compute frequency shift given theta
            if DIMWI.isFreqMW 
                freq_mw     = HCFM_obj.FrequencyMyelin(obj.x_i, obj.x_a, g, DIMWI.theta, obj.E);
                freq_mwbg   = freq_mw + fbg;
            else
                freq_mwbg   = freq_mw + fbg;
            end
            if DIMWI.isFreqIW 
                freq_iw     = HCFM_obj.FrequencyAxon(obj.x_a, g, DIMWI.theta);
                freq_iwbg   = freq_iw + fbg;
            else
                freq_iwbg   = freq_iw + fbg;
            end
            if DIMWI.isFreqEW
                freq_ew     = 0;
                freq_ewbg   = freq_ew + fbg;
            else
                freq_ewbg   = freq_ew + fbg;
            end
            %%%%%%%%%% extra decay on extracellular water estimated by HCFM %%%%%%%%%%
            if DIMWI.isR2sEW
                
                % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
                d_e = HCFM_obj.DephasingExtraaxonal(fvf, g, obj.x_i, obj.x_a, DIMWI.theta);
                
            else
                d_e = 0;
            end

            s_mw = obj.GRE_signal(Amw, -1/t2smw, freq_mwbg, te);
            s_iw = obj.GRE_signal(Aiw, -1/t2siw, freq_iwbg, te);
            s_ew = obj.GRE_signal(Aew, -1/t2sew, freq_ewbg, te) .* exp(-d_e);
    
            s = (s_mw + s_iw + s_ew) .* exp(1i*pini);

           
            % Jacobian
            if nargout == 2
                dsdAmw       = (s_mw./Amw) .* exp(1i*pini);
                dsdAiw       = (s_iw./Aiw) .* exp(1i*pini);
                dsdAew       = (s_ew./Aew) .* exp(1i*pini);
                dsdR2smw     = -te .* s_mw .* exp(1i*pini);
                dsdR2siw     = -te .* s_iw .* exp(1i*pini);
                dsdR2sew     = -te .* s_ew .* exp(1i*pini);
                dsdfreq_mwbg = 1i*2*pi*te .* s_mw .* exp(1i*pini);
                dsdfreq_iwbg = 1i*2*pi*te .* s_iw .* exp(1i*pini);
                dsdfreq_ewbg = 1i*2*pi*te .* s_ew .* exp(1i*pini);
                dsdpini      = 1i * s;

                ds = cat(1,dsdAmw(:),dsdAiw(:),dsdAew(:),dsdR2smw(:),dsdR2siw(:),dsdR2sew(:),...
                            dsdfreq_mwbg(:),dsdfreq_iwbg(:),dsdfreq_ewbg(:),dsdpini(:));
            end
        
        end
        

    end

    methods (Static)
        
        function g      = gratio(Vm, Va)

            g = sqrt( Va /(Va + Vm));

        end
        
        function s = GRE_signal(A, R2s, freq, t)
            s = A .* exp(-t .* (R2s - 1i*2*pi*freq));
        end
    end

end