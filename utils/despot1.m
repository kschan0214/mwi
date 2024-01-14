classdef despot1
% This class implements the equations in 
% Rapid combined T1 and T2 mapping using gradient recalled
% acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 4 January 2024
% Date modified:
%
    properties (Constant)
        gyro = 42.57747892;
    end

    properties (GetAccess = public, SetAccess = protected)
        tr;
        fa;

    end

    methods
        function obj = despot1(tr,fa)
            obj.tr = double(tr(:));
            obj.fa = double(fa(:));

        end

        function [t1, m0, mask_fitted] = estimate(obj,img,mask,b1)
        % [t1, m0, mask_fitted] = despot1_mapping(img,fa,tr,mask,b1)
        %
        % Input
        % --------------
        % img           : magnitude image, 4D [row,col,slice,flip_angle]
        % fa            : flip angle in degree
        % tr            : repetition time, either ms or s
        % mask          : signal mask, optional
        % b1            : B1 map, optional
        %
        % Output
        % --------------
        % t1            : T1 map, unit depends on tr
        % m0            : proton density map
        % mask_fitted   : mask of fitted voxels without encountering error
        %
        % Description: T1 mapping using DESPOT1 formalism
        %
        % Kwok-shing Chan @ DCCN
        % k.chan@donders.ru.nl
        % Date created: 30 Oct 2020
        % Date modified: 17 Nov 2020
        % Date modified: 21 July 2022
        %
            
            dims = size(img);
            
            % if B1 is not provided
            if nargin < 4 || isempty(b1)
                b1 = ones(dims(1:3));
            end
            % if mask is not provided
            if nargin < 3 || isempty(mask)
                mask = ones(dims(1:3));
            end
            
            % reshape for voxelwise operation
            img     = reshape(abs(img),[prod(dims(1:3)) dims(4)]);
            mask    = reshape(mask,[prod(dims(1:3)) 1]);
            b1      = reshape(b1,[prod(dims(1:3)) 1]);
            
            % masking
            ind = find(mask>0);
            img = img(ind,:);
            b1  = b1(ind);
            
            % DESPOT1
            t1 = zeros([prod(dims(1:3)) 1]);
            m0 = zeros([prod(dims(1:3)) 1]);
            for k = 1:size(img,1)
                if sum(img(k,:)) ~= 0
                    [m0(ind(k)), t1(ind(k))] = obj.fit(img(k,:),b1(k),'regression');
                else
                    % avoid matrix inversion with zero matrix
                    t1(ind(k)) = 0; m0(ind(k)) = 0;
                end
            end
            t1      = reshape(t1,dims(1:3));
            m0      = reshape(m0,dims(1:3));
            mask    = reshape(mask,dims(1:3));
            
            % remove physically infeasible values
            mask_fitted            = ones(size(t1));
            mask_fitted(t1<=0)     = 0;                 % T1=0 is also physically not feasible
            mask_fitted(isnan(t1)) = 0;
            mask_fitted(isinf(t1)) = 0;
            mask_fitted(m0<0)      = 0;
            mask_fitted(isnan(m0)) = 0;
            mask_fitted(isinf(m0)) = 0;
            mask_fitted            = mask_fitted .* mask;
            
            t1 = t1 .* mask_fitted;
            m0 = m0 .* mask_fitted;
 
        end

        function [m0,t1] = fit(obj,y,b1,solver,options)

            % if the fitting option is not provided then set up here
            if nargin < 5 && strcmpi(solver,'lsqnonlin')
                %%%%%%%%%% set fitting options %%%%%%%%%%
                options = optimoptions(@lsqnonlin,'MaxIter',500,...
                    'StepTolerance',1e-6,'FunctionTolerance',1e-6,'Display','off',...
                    'Jacobian','on');
            end

            % check error
            if length(y) < 2
                error('Only one point to fit.');
            end
            if length(y) ~= length(obj.fa)
                error('Size of signal is not equal to size of flip angle');
            end

            if nargin < 3 || isempty(b1)
                b1 = 1;
            end

            if nargin < 4
                solver = 'regression';
            end
            
            % Core
            % DESPOT1 formulation
            S       = double(y);
            alpha   = double(obj.fa)*b1;
            TR      = double(obj.tr);
            
            y       = S(:)./sind(alpha(:));
            xCol    = S(:)./tand(alpha(:));
            x       = ones(length(S),2);
            x(:,1)  = xCol;
                
            % solve DESPOT1 using different approach
            switch solver
                case 'regression'
                    b = x\y;
                    t1 = -TR/log(b(1));
                    m0 = b(2)/(1-exp(-TR/t1));
                case 'lsqnonneg'
                    b = lsqnonneg(x,y);
                    t1 = -TR/log(b(1));
                    m0 = b(2)/(1-exp(-TR/t1));
                case 'lsqnonlin'
                    % Obtain initial guesses
                    b = x\y;
                    t10 = -TR/log(b(1));
                    m00 = b(2)/(1-exp(-TR/t10));
                    T1_lb = 0; T1_ub = 5;
                    x0 = [t10, m00];
                    lb = [T1_lb, min(S)];
                    ub = [T1_ub, 2*m00];
            
                    res = lsqnonlin(@(x)residual(x,S,b1),x0,lb,ub,options);
                    t1 = res(2);
                    m0 = res(1);
            end
            t1 = real(t1);
            m0 = real(m0);
        end

        function [residuals, jacobian] = residual(obj, pars, y, b1)

            % forward model 
            [s,jacobian] = obj.FWD(pars, b1);

            % compute resiudal 
            residuals = s(:) - y(:);

            if ~isreal(y)
                residuals = [real(residuals);imag(residuals)];
            end

            
        end

        function [S,J] = FWD(obj,pars,b1)
            M0 = pars(1);
            T1 = pars(2);

            S  = obj.model_Bloch_1T1(M0,T1,b1);
            J  = obj.jacobian_Bloch_1T1(M0,T1,b1);
            
        end

        % Bloch non-exchanging 1-pool steady-state model
        function s = model_Bloch_1T1(obj,m0,t1,b1)
        % S = model_Bloch_1T1(obj,m0,t1,b1)
        %
        % Input:
        % ------
        %	m0: Proton density, can be arbitary unit
        % 	t1: (same unit as TR)
        %   b1: relative B1 ratio
        % Output:
        % -------
        %	S: T1-weighted signal
        %
        % Description: Calculation of signal intensity with known flip angle, T1 and TR
        %
        %   Author: Kwok-shing Chan @ University of Aberdeen
        %   Date created: Jan 1, 2016
        %   Ref: Rapid combined T1 and T2 mapping using gradient recalled
        %   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
        %   Date last edited: 20 February 2018
        %
            
            alpha = obj.fa.*b1;

            % Core algorithm
            E1 = exp(-obj.tr./t1);
            s = (m0.*(1-E1).*sind(alpha))./(1-E1.*cosd(alpha));

        end

        function j = jacobian_Bloch_1T1(obj,m0,t1,b1)

            alpha = obj.fa.*b1;

            % Core algorithm
            E1 = exp(-obj.tr./t1);

            dSdM0 = (1-E1).*sind(alpha)./(1-E1.*cosd(alpha));
            dSdT1 = m0.*sind(alpha).*(cosd(alpha)-1)*exp(obj.tr/t1).*obj.tr./(t1.^2 .* (cosd(alpha)-exp(obj.tr./t1)).^2);
            j = [dSdM0(:), dSdT1(:)];
        end
   
    end

    methods(Static)

        % Calculation of Ernst angle with knwon T1 and TR
        function EA = ernst_angle(T1,TR)
        % function EA = ernst_angle(T1,TR)
        %   Input:
        %           -   T1: (same unit as TR)
        %           -   TR: (same unit as T1)
        %   Output:
        %           -   EA: Ernst angle based on T1 and TR
        %
        %   Author: Kwok-shing Chan @ University of Aberdeen
        %   Date created: Jan 1, 2016
        %   Ref: Rapid combined T1 and T2 mapping using gradient recalled
        %   acquisition in the steady state, Deoni et al. MRM 2003;49:515-526
        %
            EA = acosd(exp(-TR/T1));
        end
    end

end