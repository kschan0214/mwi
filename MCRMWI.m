classdef MCRMWI
    properties
        gyro    = 42.57747892;
        B0      = 3; % T
        B0dir   = [0;0;1]; % main magnetic field direction with respect to FOV
        % tissue parameters
        x_i     = -0.1;     % ppm
        x_a     = -0.1;     % ppm
        E       = 0.02;     % ppm
        rho_mw  = 0.42;     % relative water proton density relative to IEW
        t1_mw   = 234e-3;   % s
        te;
        tr;
        fa;
        epgx_params;

    end

    methods
        %% Constructor
        function obj = MCRMWI(te,tr,fa,epgx_params,fixed_params)
            obj.te = double(te(:));
            obj.tr = double(tr(:));
            obj.fa = double(fa(:));

            epgx_params_default.isExchange = false;
            epgx_params_default.isEPG      = false;
            epgx_params_default.npulse     = 200;
            epgx_params_default.rfphase    = 50;
            
            % EPG-X options
            if nargin >=4
                if isfield(epgx_params,'isExchange')
                    epgx_params_default.isExchange = epgx_params.isExchange;
                end
                if isfield(epgx_params,'isEPG') 
                    epgx_params_default.isEPG = epgx_params.isEPG;
                end
                if isfield(epgx_params,'npulse')
                    epgx_params_default.npulse = epgx_params.npulse;
                end
                if isfield(epgx_params,'rfphase')
                    epgx_params_default.rfphase = epgx_params.rfphase;
                end
            end
            obj.epgx_params    = epgx_params_default;
            
            % fixed tissue and scanner parameters
            if nargin == 5 
                if isfield(fixed_params,'x_i')
                    obj.x_i     = double(fixed_params.x_i);
                end
                if isfield(fixed_params,'x_a')
                    obj.x_a     = double(fixed_params.x_a);
                end
                if isfield(fixed_params,'E')
                    obj.E       = double(fixed_params.E);
                end
                if isfield(fixed_params,'rho_mw')
                    obj.rho_mw  = double(fixed_params.rho_mw);
                end
                if isfield(fixed_params,'B0')
                    obj.B0      = double(fixed_params.B0);
                end
                if isfield(fixed_params,'B0dir')
                    obj.B0dir   = double(fixed_params.B0dir);
                end
                if isfield(fixed_params,'t1_mw')
                    obj.t1_mw   = double(fixed_params.t1_mw);
                end
            end

        end

        %% Data fitting
        % prepare data for MCR-MWI fitting
        function [algoPara,data_obj_all] = data_preparation(obj,algoPara,imgPara)

            disp('=================================================');
            disp('Myelin water imaing: MCR-(DI)MWI Data Preparation');
            disp('=================================================');

            %%%%%%%%%% create directory for temporary results %%%%%%%%%%
            default_output_filename = 'mcrmwi_results.mat';
            temp_prefix             = 'temp_mcrmwi_';
            [output_dir, output_filename]       = obj.setup_output(imgPara,default_output_filename);
            [temp_dir,temp_prefix, identifier]  = obj.setup_temp_dir(imgPara,temp_prefix);

            %%%%%%%%%% validate algorithm and image parameters %%%%%%%%%%
            [algoPara,imgPara] = obj.check_and_set_default(algoPara,imgPara);

            %%%%%%%%%% capture all image parameters %%%%%%%%%%
            data  = double(imgPara.img);
            b1map = double(imgPara.b1map);
            mask  = imgPara.mask>0;
            fm0   = double(imgPara.fieldmap);
            pini0 = double(imgPara.pini);

            %%%%%%%%%% DIMWI %%%%%%%%%%
            [icvf, ff, theta] = obj.setup_DIMWI(imgPara, algoPara);

            %%%%%%%%%% prepare input %%%%%%%%%%
            % parfor could be slowed down by a large number of unmasked voxels, compute a new mask based on input data
            mask_valid_DIMWI = squeeze(sum(ff,4))>0;
            % mask    = and(mask>0,squeeze(sum(ff,4))>0);
            if ~algoPara.DIMWI.isFitVic
                mask_valid_DIMWI = and(mask_valid_DIMWI,icvf>0);
                % mask = and(mask>0,icvf>0);
            end
            
            if algoPara.isNormData
                [scaleFactor, data] = mwi_image_normalisation(data, mask);
            else
                scaleFactor = 1;
            end

            %%%%%%%%%% check advanced starting point strategy %%%%%%%%%%
            % use lowest flip angle to estimate R2*
            [~,ind]         = min(obj.fa);
            r2s0            = R2star_trapezoidal(data(:,:,:,:,ind),obj.te);
            % r2s0            = R2star_trapezoidal(mean(data(:,:,:,:,:),5),obj.te);
            mask_valida_r2s = and(~isnan(r2s0),~isinf(r2s0));
            r2s0(mask_valida_r2s == 0) = 0;

            % only works for 3T data
            advancedStarting = algoPara.advancedStarting;
            if strcmpi(advancedStarting,'default') || strcmpi(advancedStarting,'robust')
                fprintf('Estimate starting points using predefined model...')
                
                if obj.B0 > 2.5 && obj.B0 < 3.5 % 3T
                    t2s_pre = [10e-3,60e-3];    % [T2sMW, T2sIEW] in second
                    t1_pre  = [234e-3, 1];    	% [T1MW, IEW], in second
                elseif obj.B0 > 1 && obj.B0 < 2 % 1.5T
                    t2s_pre = [10e-3,80e-3];    
                    t1_pre  = [234e-3, 0.7];    
                elseif obj.B0 > 6 && obj.B0 < 2 % 7T
                    t2s_pre = [10e-3,40e-3];    
                    t1_pre  = [234e-3, 1.5];    
                end
                
                switch advancedStarting
                    case 'default'
                        [m00,mwf0,t2siew0,t1iew0] = superfast_mwi_2m_mcr_self(data,obj.te,obj.fa,obj.tr,t2s_pre,t1_pre(1),mask,b1map,'superfast');
                    case 'robust'
                        [m00,mwf0] = superfast_mwi_2m_mcr(data,obj.te,obj.fa,obj.tr,t2s_pre,t1_pre,mask,b1map);
                end
                m00 = sum(m00,4); % total water
                % also masked out problematic voxels detected by superfast method
                mask_valid_m0 = m00 > 0;
                disp('Completed.')
            end
            % Simple DESPOT1 T1 mapping
            if ~exist('t1iew0','var')
                fprintf('Estimate T1 using DESPOT1...')
                [t1iew0, m0] = despot1_mapping(permute(data(:,:,:,1,:),[1 2 3 5 4]),obj.fa,obj.tr,mask,b1map);
                if ~exist('m00','var')
                    m00 = m0;
                    % also masked out DESPOT1 problematic voxels
                    mask_valid_m0 = m00 > 0;
                end
                disp('Completed.')
            end
            
            % final mask for fitting
            mask = and(and(and(mask,mask_valida_r2s),mask_valid_m0),mask_valid_DIMWI);
            
            %%%%%%%%%% Batch processing preparation %%%%%%%%%%
            % set up data_obj for batch processing
            % input data
            data_obj_all = setup_batch_create_data_obj_slice(data,      mask, 'data');
            data_obj_all = setup_batch_create_data_obj_slice(b1map,     mask, 'b1map',  data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(fm0,       mask, 'fm0',    data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(pini0,     mask, 'pini0',  data_obj_all);
            % DIMWI
            data_obj_all = setup_batch_create_data_obj_slice(icvf,      mask, 'icvf',   data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(theta,     mask, 'theta',  data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(ff,        mask, 'ff',     data_obj_all);
            % derived initial guess
            data_obj_all = setup_batch_create_data_obj_slice(m00,       mask, 'm00',	data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(t1iew0,    mask, 't1iew0',	data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(r2s0,      mask, 'r2s0',	data_obj_all);
            if exist('mwf0','var');     data_obj_all = setup_batch_create_data_obj_slice(mwf0,      mask, 'mwf0',       data_obj_all); end
            if exist('t2siew0','var');  data_obj_all = setup_batch_create_data_obj_slice(t2siew0,   mask, 't2siew0',	data_obj_all); end
            
            for kbat = 1:numel(data_obj_all)
                data_obj    = data_obj_all(kbat);
                batchNumber = kbat;
                save(fullfile(temp_dir,strcat(temp_prefix,'_batch',num2str(kbat))),'data_obj','identifier','scaleFactor','batchNumber');
            end
            algoPara.identifier         = identifier;
            algoPara.temp_dir           = temp_dir;
            algoPara.temp_prefix        = temp_prefix;
            algoPara.output_dir         = output_dir;
            algoPara.output_filename    = output_filename;

            save(output_filename,'algoPara');

            disp('Note that the final mask could be different from the input mask due to data validity.');
            disp('Data preparation step is completed.');

        end

        % Data fitting on whole dataset
        function res_obj = estimate(obj,algoPara,data_obj)

            disp('==============================================');
            disp('Myelin water imaing: MCR-(DI)MWI model fitting');
            disp('==============================================');

            % display some messages
            obj.display_message(algoPara);
    
            % determine if the data is provided directly or saved to the disk
            if nargin < 3
                % if no data input then get identifier from algoPara
                data_obj    = [];
                temp_files  = fullfile(algoPara.temp_dir,algoPara.temp_prefix);
                filelist    = dir([temp_files '*']);
                if ~isempty(filelist)
                    Nbatch = numel(filelist);
                end

                isBatch = true;

            elseif ~isempty(data_obj)
                isBatch = false;
                Nbatch  = 1;
            end

            %%%%%%%%%% create directory for temporary results %%%%%%%%%%
            default_output_filename = 'mwi_mcr_results.mat';
            [output_dir, output_filename]                = obj.setup_output(algoPara,default_output_filename);
            identifier = algoPara.identifier;

            %%%%%%%%%% log command window display to a text file %%%%%%%%%%
            logFilename = fullfile(output_dir, ['run_mwi_' identifier '.log']);
            logFilename = check_unique_filename_seqeunce(logFilename);
            diary(logFilename)

            fprintf('Output directory                : %s\n',output_dir);
            fprintf('Intermediate results identifier : %s\n',identifier);

            %%%%%%%%%% capture all algorithm parameters %%%%%%%%%%
            maxIter         = algoPara.maxIter;
            fcnTol          = algoPara.fcnTol;
            stepTol         = algoPara.stepTol;
            isParallel      = algoPara.isParallel;
            numEst          = algoPara.numEst;

            % data fitting method related parameters
            fitParams.isComplex      = algoPara.isComplex;
            fitParams.isWeighted     = algoPara.isWeighted;
            fitParams.weightMethod   = algoPara.weightMethod;
            fitParams.weightPower    = algoPara.weightPower;
            fitParams.isFastEPG      = algoPara.isFastEPG;
            fitParams.userDefine     = algoPara.userDefine;
            fitParams.isInvivo       = algoPara.isInvivo;
            fitParams.isSelfBoundary = algoPara.isSelfBoundary;
            fitParams.isFitT1mw      = algoPara.isFitT1mw;
            fitParams.DEBUG          = algoPara.DEBUG;

            % DIMWI parameters
            DIMWI = algoPara.DIMWI;

            %%%%%%%%%% set fitting options %%%%%%%%%%
            options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*numEst,...
                'StepTolerance',stepTol,'FunctionTolerance',fcnTol,'Display','off');

            if fitParams.DEBUG
                % if DEBUG is on then disables parallel computing
                isParallel      = false;
                options.Display = 'iter';
            end

            % create a parallel pool
            if isParallel
                gcp;
            end

            % main
            try
                
                % %%%%%%%%%% initiate progress display %%%%%%%%%%
                % fprintf('%i voxel(s) to be fitted...\n',numMaskedVoxel);
                % if exist(temp_filename,'file')
                %     % restore progress
                %     disp('Restoring previous progress...')
                %     load(temp_filename);
                %     isRestore = true;
                % else
                %     % new progress
                %     fbat = 0;
                %     lastBatchEndTime = nElement *8; %roughly 8s for 1 voxel
                %     isRestore = false;
                % end
                % progress_display(numBatch,fbat,lastBatchEndTime,isRestore);
                % isRestore = false;

                %%%%%%%%%% fitting main %%%%%%%%%%
                for kbat = 1:Nbatch
                    if isBatch
                        % if batch mode then check if the input mat
                        % contains output variable
                        temp_mat = strcat(temp_files,['_batch' num2str(kbat)]);
                        variableInfo = who('-file', temp_mat);
                        isFit = ~ismember('res_obj', variableInfo);
                    else
                        isFit = true;
                    end
                    
                    if isFit 
                        if isBatch
                            % if batch mode then load data from disk
                            load(temp_mat);
                        end

                        % process non-empty slices
                        if ~data_obj.isEmptySlice

                            data    = double(data_obj.data);
                            b1map   = double(data_obj.b1map);
                            
                            extra_data = [];
                            for k = 1:size(data,1)

                                % DIMWI
                                extra_data(k).icvf  = double(data_obj.icvf(k,:));      % Intra-neurite fraction
                                extra_data(k).ff    = double(data_obj.ff(k,:));        % fibre fraction
                                extra_data(k).theta = double(data_obj.theta(k,:));     % angle between B0 and fibre

                                % initial starting points
                                extra_data(k).pini0 = double(data_obj.pini0(k,:)); % initial phase
                                extra_data(k).db0   = double(data_obj.fm0(k,:));   % fieldmap
                                extra_data(k).m00   = double(data_obj.m00(k));   % proton density weighted isgnal
                                   
                                extra_data(k).r2s0  = double(data_obj.r2s0(k));   % R2*
                                extra_data(k).t1iew0= double(data_obj.t1iew0(k)); 
                                
                                if isfield(data_obj,'mwf0');    extra_data(k).mwf0      = double(data_obj.mwf0(k));         end % MWF
                                if isfield(data_obj,'t2siew0'); extra_data(k).t2siew0   = double(data_obj.t2siew0(k));      end
                                if isfield(data_obj,'t1iew0');  extra_data(k).t1iew0    = double(data_obj.t1iew0(k));       end % T1 IEW
                            end
                            
                            % start timer
                            tic;
                            % create an empty array for fitting results 
                            resnorm   = zeros(size(data,1),1);
                            iter      = zeros(size(data,1),1);
                            exitflag  = zeros(size(data,1),1);
                            if isParallel
                                parfor k = 1:size(data,1)
                                    s           = squeeze(data(k,:,:));
                                    b1          = b1map(k);
                                    extraData   = extra_data(k);
                                    
                                        [fitRes(k), resnorm(k),exitflag(k),output] = ...
                                            obj.fit(s,b1,extraData,DIMWI,fitParams,options);
                                        iter(k) = output.iterations;
                                end
                            else
                                for k = 1:size(data,1)
                                    s           = squeeze(data(k,:,:));
                                    b1          = b1map(k);
                                    extraData   = extra_data(k);
                                    
                                        [fitRes(k), resnorm(k),exitflag(k),output] = ...
                                            obj.fit(s,b1,extraData,DIMWI,fitParams,options);
                                        iter(k) = output.iterations;
                                end
                            end
                            lastBatchEndTime = toc;

                            % reshape the fitting results into image
                            ind = find(data_obj.mask(:)>0);
                            fieldname = fieldnames(fitRes);
                            for kfield = 1:numel(fieldname)
                                res_obj.(fieldname{kfield}) = zeros([numel(data_obj.mask(:))  numel(fitRes(1).(fieldname{kfield}))]);
                            end
                            for kfield = 1:numel(fieldname)
                                for k = 1:numel(ind)
                                    res_obj.(fieldname{kfield})(ind(k),:) = fitRes(k).(fieldname{kfield});
                                end   
                            end
                            for kfield = 1:numel(fieldname)
                                res_obj.(fieldname{kfield}) = reshape(res_obj.(fieldname{kfield}),[size(data_obj.mask),size(res_obj.(fieldname{kfield}),2)]);
                            end
                            
                            rn = zeros(size(data_obj.mask)); rn(ind) = resnorm;
                            it = zeros(size(data_obj.mask)); it(ind) = iter;
                            ef = zeros(size(data_obj.mask)); ef(ind) = exitflag;
                            
                            res_obj.resnorm      = rn;
                            res_obj.iterations   = it;
                            res_obj.exitflag     = ef;
                            
                            % display progress
                            % progress_display(numBatch,fbat,lastBatchEndTime,isRestore);
                            if isBatch
                                save(temp_mat,'res_obj','-append')
                            end
                        end
                    end
                end

            catch ME
    
                % close log file
                disp('There was an error! Please check the command window/error message file for more information.');
                diary off
                
                % open a new text file for error message
                errorMessageFilename = fullfile(output_dir, ['run_mwi_' identifier '.error']);
                errorMessageFilename = check_unique_filename_seqeunce(errorMessageFilename);
                fid = fopen(errorMessageFilename,'w');
                fprintf(fid,'The identifier was:\n%s\n\n',ME.identifier);
                fprintf(fid,'The message was:\n\n');
                msgString = getReport(ME,'extended','hyperlinks','off');
                fprintf(fid,'%s',msgString);
                fclose(fid);
                
                % rethrow the error message to command window
                rethrow(ME);
            end

        end

        % Data fitting on 1 voxel
        function [fitRes, resnorm,exitflag,output] = fit(obj,y,b1,extraData, DIMWI, fitParams,options)

            % if the fitting option is not provided then set up here
            if nargin < 7
                %%%%%%%%%% set fitting options %%%%%%%%%%
                options = optimoptions(@lsqnonlin,'MaxIter',fitParams.maxIter,'MaxFunctionEvaluations',200*fitParams.numEst,...
                    'StepTolerance',fitParams.stepTol,'FunctionTolerance',fitParams.fcnTol,'Display','off');
            end

            %%%%%%%%%% Step 1: Prepare fitting algorithm %%%%%%%%%%
            % set the starting point, lower bound and upper bound
            [x0,lb,ub] = obj.set_up_x0lbub(y,extraData,DIMWI,fitParams);

            % precompute EPG-X's transition matrix here for speed
            Nfa     = numel(obj.fa);
            T3D_all = cell(1,Nfa);
            if obj.epgx_params.isExchange && obj.epgx_params.isEPG
                phiCycle = RF_phase_cycle(obj.epgx_params.npulse,obj.epgx_params.rfphase);
                for kfa=1:Nfa
                    T3D_all{kfa} = obj.PrecomputeT(phiCycle,d2r(obj.fa(kfa)*b1));
                end
            end

            %%%%%%%%%% Step 2: compute the weighting terms %%%%%%%%%%
            if fitParams.isWeighted
                switch lower(fitParams.weightMethod)
                    case 'norm'
                        % weights using echo intensity, as suggested in Nam's paper
                        w = sqrt(abs(y));
                    case '1stecho'
                        p = fitParams.weightPower;
                        % weights using the 1st echo intensity of each flip angle
                        w = bsxfun(@rdivide,abs(y).^p,abs(y(1,:)).^p);
                end
            else
                % compute the cost without weights
                w = ones(size(y));
            end

            %%%%%%%%%% Step 3: run fitting algorithm  %%%%%%%%%%
            DIMWI.theta = extraData.theta;
            DIMWI.ff    = extraData.ff;
            DIMWI.icvf  = extraData.icvf;
            [x,resnorm,~,exitflag,output] = lsqnonlin(@(x)obj.residuals(x,y,b1,DIMWI, w,fitParams,T3D_all),x0,lb,ub,options);

            %%%%%%%%%% Step 4: get fitted parameters  %%%%%%%%%%
            counter = 1;
            fitRes.S0_MW = x(counter);          counter = counter +1;
            if ~DIMWI.isFitVic
                fitRes.S0_IEW = x(counter);     counter = counter +1;
            else
                fitRes.S0_IW = x(counter);      counter = counter +1;
                fitRes.S0_EW = x(counter);      counter = counter +1;
            end
            fitRes.R2s_MW = x(counter);         counter= counter + 1;
            fitRes.R2s_IW = x(counter);         counter= counter + 1;
            if DIMWI.isFitR2sEW
                fitRes.R2s_EW = x(counter);     counter= counter + 1;
            end
            if fitParams.isFitT1mw
                fitRes.T1_MW = x(counter);      counter= counter + 1;
            end
            fitRes.T1_IEW = x(counter);         counter= counter + 1;
            if obj.epgx_params.isExchange
                fitRes.kiewm = x(counter);      counter= counter + 1;
            end
            if DIMWI.isFitFreqMW
                fitRes.Freq_MW = x(counter);    counter= counter + 1;
            end
            if DIMWI.isFitFreqIW
                fitRes.Freq_IW = x(counter);    counter= counter + 1;
            end
            if fitParams.isComplex
                fitRes.Freq_BKG = x(counter:counter+length(obj.fa)-1); 
                fitRes.pini = x(counter+length(obj.fa):end);
            end

        end

        % compute the residuals and jacobian for data fitting
        function [residuals, jacobian] = residuals(obj, pars, y, b1,  DIMWI, w,fitParams, T3D_all)

            s = obj.FWD(pars, b1, T3D_all, DIMWI, fitParams);

            residuals = s(:) - y(:);
            residuals = residuals .* w(:);

            if fitParams.isComplex
                residuals = [real(residuals);imag(residuals)];
            end
            
        end

        % MCR-(DI)MWI Forward model fot data fitting
        function s = FWD(obj, pars, b1, T3D_all, DIMWI, fitParams)

            EPGXParams = obj.epgx_params;

            % number of flip angles and echo times
            Nfa = numel(obj.fa);
            Nte = numel(obj.te);

            Amw = pars(1);

            if ~DIMWI.isFitVic
                Aie = pars(2); 
                Aiw = Aie * DIMWI.icvf; 
                Aew = Aie * (1-DIMWI.icvf);
                counter = 3;
            else
                Aiw = pars(2); 
                Aew = pars(3); 
                counter = 4;
            end

            % T2*s
            t2smw = 1/pars(counter);        counter = counter + 1;
            t2siw = 1/pars(counter);        counter = counter + 1;
            if DIMWI.isFitR2sEW
                t2sew = 1/pars(counter);    counter = counter + 1;
            else
                t2sew = t2siw;
            end
            % T1
            if fitParams.isFitT1mw
                t1mw = pars(counter);       counter = counter + 1;
            else
                t1mw = obj.t1_mw; 
            end
            t1iew = pars(counter);          counter = counter + 1;
            if EPGXParams.isExchange
                kiewmw = pars(counter);     counter = counter + 1;
            else
                kiewmw = 0;
            end
            % frequency shifts
            if DIMWI.isFitFreqMW
                freq_mw = pars(counter); counter = counter + 1;
            else
                freq_mw = 0;
            end
            if DIMWI.isFitFreqIW
                freq_iw = pars(counter); counter = counter + 1;
            else
                freq_iw = 0;
            end

            % external effects
            if ~fitParams.isComplex % magnitude fitting
                fbg     = zeros(Nfa);                          
                pini    = 0;
            else    % other fittings
                fbg     = pars(counter:counter+Nfa-1);  
                pini    = pars(counter+Nfa:end);
            end
            
            freq_ew = 0;

            %%%%%%%%%% simulate signal based on parameter input %%%%%%%%%%
            s = zeros([Nte Nfa length(DIMWI.ff)]);
            DIMWI_curr = DIMWI;
            for kfo = 1:length(DIMWI.theta)
                
                s(:,:,kfo) = obj.model_MCRMWI(b1,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1mw,t1iew,kiewmw,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI_curr,T3D_all);
            end
            s = sum(bsxfun(@times,s,permute(DIMWI.ff(:),[2 3 1])),3);
            
        end

        %% Signal models
        % Bloch-McConnell 2-pool exchange steady-state signal
        function [s,ss_pool] = model_BM_2T1(obj,M0r,M0f,t1r,t1f,krf,kfr,b1)
            if nargin < 8
                b1 = 1;
            end

            % B1 corrected flip angle
            alpha = obj.fa * b1;

            % no. of flip angles
            Nfa = length(obj.fa);
            
            % derive the rest of the exchange rates
            if isempty(krf)
                if M0r ~= 0
                    krf = kfr * M0f/M0r;
                else
                    krf = 0;
                end
                    
            else
                if M0f ~= 0
                    kfr   = krf * M0r/M0f;
                else
                    kfr = 0;
                end
            end
            
            % convert T1s to R1s for simplicity
            R1r    = 1/t1r;
            R1f    = 1/t1f;
            
            %%% Relaxation-exchange matrix for Longitudinal components
            Lambda_L = [-R1r-krf    , kfr       ;
                        krf      	, -R1f-kfr  ];
                    
            Xi_L = expm(obj.tr*Lambda_L);
            
            % Use same notation as paper
            C_L = [R1r*M0r; R1f*M0f];
            I = eye(length(C_L));
                
            ss_pool = zeros(length(C_L),Nfa);
            s = zeros(1,Nfa);
            for kfa = 1:Nfa
                
                % Use same notation as paper
                T = cosd(alpha(kfa))*I;
                
                % Spencer and Fishbein J. Magn. Reson:142(2000)120-135
                ss_pool(:,kfa) = sind(alpha(kfa))* ((I-Xi_L*T) \ (Xi_L - I)) * (Lambda_L \ C_L);
                
                s(kfa) = ones(1,length(C_L)) * ss_pool(:,kfa);
            
            end
        
        end

        % Bloch non-exchanging 3-pool steady-state model
        function s = model_Bloch_3T1(obj,m0mw,m0iw,m0ew,t1mw,t1iw,t1ew,b1)

            if nargin < 8
                b1 = 1;
            end
            
            s = zeros(3,length(obj.fa));
            
            alpha = obj.fa * b1;
            
            s(1,:) = m0mw .* sind(alpha) .* (1-exp(-obj.tr/t1mw))./(1-cosd(alpha).*exp(-obj.tr/t1mw));
            s(2,:) = m0iw .* sind(alpha) .* (1-exp(-obj.tr/t1iw))./(1-cosd(alpha).*exp(-obj.tr/t1iw));
            s(3,:) = m0ew .* sind(alpha) .* (1-exp(-obj.tr/t1ew))./(1-cosd(alpha).*exp(-obj.tr/t1ew));
    
    
        end

        % Bloch non-exchanging 1-pool steady-state model
        function s = model_Bloch_1T1(obj,m0,t1,b1)
            
            alpha = obj.fa.*b1;

            % Core algorithm
            E1 = exp(-obj.tr./t1);
            s = (m0.*(1-E1).*sind(alpha))./(1-E1.*cosd(alpha));

        end

        % MCR-(DI)MWI T1-T2* signal model
        function s = model_MCRMWI(obj,b1,Amw,Aiw,Aew,t2smw,t2siw,t2sew,t1mw,t1iew,kiewmw,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI,T3D_all)
            
            % number of flip angles
            Nfa = numel(obj.fa);

            % Step 1: validate input variables
            EPGX        = obj.epgx_params;

            if ~DIMWI.isFitFreqMW || ~DIMWI.isFitFreqIW || ~DIMWI.isFitR2sEW
                hcfmObj = HCFM(obj.te,obj.B0);
            end
            
            % if no phase offset set then no phase offset
            if nargin < 16
                pini=0;
            end
            
            % replicate phase offset for all acquisition if only 1 value is available
            if numel(pini) == 1
                pini = pini * ones(Nfa,1);
            end
            
            % replicate total field for all acquisition if only 1 value is available
            if numel(fbg) == 1
                fbg = fbg * ones(Nfa,1);
            end
            
            % Step 2: get DIMWI derived parameters
            %%%%%%%%%% frequency shifts estimated using HCFM %%%%%%%%%%
            if ~DIMWI.isFitFreqMW || ~DIMWI.isFitFreqIW
                
                % derive g-ratio 
                g = hcfmObj.gratio(abs(Aiw),abs(Amw)/obj.rho_mw);
                
                % compute frequency shifts given theta
                if ~DIMWI.isFitFreqMW 
                    freq_mw = hcfmObj.FrequencyMyelin(obj.x_i,obj.x_a,g,DIMWI.theta,obj.E);
                end
                if ~DIMWI.isFitFreqMW 
                    freq_iw = hcfmObj.FrequencyAxon(obj.x_a,g,DIMWI.theta);
                end
            end
            
            %%%%%%%%%% extra decay on extracellular water estimated by HCFM %%%%%%%%%%
            if ~DIMWI.isFitR2sEW
                
                % assume extracellular water has the same T2* as intra-axonal water
                t2sew = t2siw;

                fvf = hcfmobj.FibreVolumeFraction(abs(Aiw),abs(Aiw),abs(Amw)/obj.rho_mw);
                % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
                d_e = hcfmobj.DephasingExtraaxonal(fvf,g,obj.x_i,obj.x_a,DIMWI.theta);
                
                
            else
                d_e = zeros(size(obj.te));
            end
            
            % Compute steady-state signal for each pool
            % 1: MW; 2: IW; 3: EW
            if EPGX.isExchange
                % exchange assumed to take place between (myelin water + myelin macromolecules) and IEW, except when EPGX.rho_mw = 1
                
                % derive tissue properties
                myelinVolumeSignal      = Amw/obj.rho_mw;
                IEWVolumeSignal         = Aiw + Aew;
                volumeScaleFactor       = IEWVolumeSignal + myelinVolumeSignal;
                v_ic                    = Aiw / IEWVolumeSignal;
                    
                if EPGX.isEPG
                    
                    % EPG-X
                    phiCycle    = RF_phase_cycle(EPGX.npulse,EPGX.rfphase);
                    t1x         = [t1iew, t1mw]; 
                    t2x         = [t2siw,t2smw]; % assuming T2* of iw has similar T2 of long T1 compartment
            
                    fx = myelinVolumeSignal/volumeScaleFactor;  % myelin volume fraction
                    fs = (freq_mw-freq_iw);                     % frequency difference between long and short T1 compartments
                    
                    % compute saturation factor
                    SF = zeros(length(obj.fa),2);
                    % start with steady-state signals, longitudinal magnetisation (Mz)
                    z1_all = obj.model_Bloch_1T1(1-fx,  t1iew, b1)./sind(obj.fa*b1);
                    z2_all = obj.model_Bloch_1T1(fx,    t1mw,  b1)./sind(obj.fa*b1);
                    for ii=1:length(obj.fa)

                        % true flip angle
                        alpha = obj.fa(ii)*b1;
            
                        % 2 pools, with exchange
                        z1 = z1_all(ii);
                        z2 = z2_all(ii);
                        
                        % EPG-X core
                        tmp = EPGX_GRE_BMsplit_PrecomputedT(T3D_all{ii},phiCycle,obj.tr,t1x,t2x,fx,kiewmw,'delta',fs,'kmax',10,'ss',[z1,z2]);
            
                        % saturation factors, 1: IEW; 2:Myelin  
                        SF(ii,1) = tmp{2}(end);
                        SF(ii,2) = tmp{1}(end); 
                        % SF(ii,1) = abs(tmp{2}(end)); % updated: 20231223
                        % SF(ii,2) = abs(tmp{1}(end)); 
                    end
                    ss(1,:) = SF(:,1) * volumeScaleFactor * obj.rho_mw;        % MW = Myelin volume * MW density
                    ss(2,:) = SF(:,2) * volumeScaleFactor * v_ic;               % IW = IEW * v_ic
                    ss(3,:) = SF(:,2) * volumeScaleFactor * (1 - v_ic);         % EW = IEW * (1-v_ic)
                    
                else

                    % Exchange only, no need to multiply volumeScaleFactor as direct
                    % volume is used instead of volume fraction
                    % 20220209: bug fix 
                    [~,ss_pool] = obj.model_BM_2T1(myelinVolumeSignal,IEWVolumeSignal,t1mw,t1iew,[],kiewmw,b1);
                    ss(1,:) = ss_pool(1,:) * obj.rho_mw;	     % MW = Myelin volume * MW density
                    ss(2,:) = ss_pool(2,:) * v_ic;           % IW = IEW * v_ic
                    ss(3,:) = ss_pool(2,:) * (1-v_ic);       % EW = IEW * (1-v_ic)
                     
                end
            else
                % non-exchange steady-state
                ss = obj.model_Bloch_3T1(Amw,Aiw,Aew,t1mw,t1iew,t1iew,b1);

            end
            
            % create 2D matrices for direct multiplication
            [te2D,Amw2D]    = ndgrid(obj.te,ss(1,:));
            [~,Aiw2D]       = ndgrid(obj.te,ss(2,:));
            [~,Aew2D]       = ndgrid(obj.te,ss(3,:));
            [~,fbg2D]       = ndgrid(obj.te,fbg);
            [~,pini2D]      = ndgrid(obj.te,pini);
            [d_e2D,~]       = ndgrid(d_e,obj.fa);
            
            %% T2* decay effect
            
            s = (Amw2D .* exp(te2D * (-1/t2smw + 1i*2*pi*freq_mw)) + ...                % S_MW
                 Aiw2D .* exp(te2D * (-1/t2siw + 1i*2*pi*freq_iw)) + ...                % S_IW
                 Aew2D .* exp(te2D * (-1/t2sew + 1i*2*pi*freq_ew)).*exp(-d_e2D)) .* ... % S_EW
                 exp(1i*pini2D) .* exp(1i*2*pi*fbg2D.*te2D);                            % phase offset and total field
            
        end

        %% Utility functions
        % check and set default
        function [algoPara2,imgPara2] = check_and_set_default(obj,algoPara,imgPara)
        
            imgPara2    = imgPara;
            algoPara2   = algoPara;
            
            %%%%%%%%%% 1. check algorithm parameters %%%%%%%%%%
            % check debug
            try algoPara2.DEBUG             = algoPara.DEBUG;         	catch; algoPara2.DEBUG = false; end
            % check parallel computing 
            try algoPara2.isParallel        = algoPara.isParallel;   	catch; algoPara2.isParallel = false; end
            % fast fitting when EPG enabled
            try algoPara2.isFastEPG       	= algoPara.isFastEPG;    	catch; algoPara2.isFastEPG = false; end
            % check maximum iterations allowed
            try algoPara2.maxIter           = algoPara.maxIter;     	catch; algoPara2.maxIter = 200; end
            % check function tolerance
            try algoPara2.fcnTol            = algoPara.fcnTol;      	catch; algoPara2.fcnTol = 1e-5; end
            % check step tolerance
            try algoPara2.stepTol           = algoPara.stepTol;     	catch; algoPara2.stepTol = 1e-5; end
            % check normalised data before fitting
            try algoPara2.isNormData        = algoPara.isNormData;  	catch; algoPara2.isNormData = true; end
            % check weighted sum of cost function
            try algoPara2.isWeighted        = algoPara.isWeighted;  	catch; algoPara2.isWeighted = true; end
            % check weighted sum of cost function
            try algoPara2.weightPower       = algoPara.weightPower;  	catch; algoPara2.weightPower = 2; end
            % check # of phase-corrupted echoes
            try algoPara2.isComplex         = algoPara.isComplex;       catch; algoPara2.isComplex = false; end
            % check # of phase-corrupted echoes
            try algoPara2.isInvivo          = algoPara.isInvivo;       	catch; algoPara2.isInvivo = true; end
            % check # of phase-corrupted echoes
            try algoPara2.isSelfBoundary 	= algoPara.isSelfBoundary;	catch; algoPara2.isSelfBoundary = false; end
    
            % check user bounds and initial guesses
            try algoPara2.userDefine.x0     = algoPara.userDefine.x0;   catch; algoPara2.userDefine.x0 = [];end
            try algoPara2.userDefine.lb     = algoPara.userDefine.lb;   catch; algoPara2.userDefine.lb = [];end
            try algoPara2.userDefine.ub     = algoPara.userDefine.ub;   catch; algoPara2.userDefine.ub = [];end
            % check hollow cylinder fibre model parameters
            try algoPara2.DIMWI.isFitFreqMW	= algoPara.DIMWI.isFitFreqMW;	catch; algoPara2.DIMWI.isFitFreqMW = false;end
            try algoPara2.DIMWI.isFitFreqIW	= algoPara.DIMWI.isFitFreqIW;	catch; algoPara2.DIMWI.isFitFreqIW = false;end
            try algoPara2.DIMWI.isFitR2sEW  = algoPara.DIMWI.isFitR2sEW;	catch; algoPara2.DIMWI.isFitR2sEW  = false;end
            try algoPara2.DIMWI.isFitVic    = algoPara.DIMWI.isFitVic;      catch; algoPara2.DIMWI.isFitVic    = false;end
           
            % check advanced starting points strategy
            try algoPara2.advancedStarting  = algoPara.advancedStarting;catch; algoPara2.advancedStarting = [];end
            
            %%%%%%%%%% 2. check data integrity %%%%%%%%%%
            disp('-----------------------');
            disp('Checking data integrity');
            disp('-----------------------');
            % check if the number of echo times matches with the data
            if numel(obj.te) ~= size(imgPara.img,4)
                error('The length of TE does not match with the 4th dimension of the image.');
            end
            if numel(obj.fa) ~= size(imgPara.img,5)
                error('The length of flip angle does not match with the 5th dimension of the image.');
            end
            % check signal mask
            try
                imgPara2.mask = imgPara.mask;
                disp('Mask input                : True');
            catch
                imgPara2.mask = max(max(abs(imgPara.img),[],4),[],5)./max(abs(imgPara.img(:))) > 0.05;
                disp('Mask input                : false');
                disp('Default masking method is used.');
            end
            % check field map
            try
                imgPara2.fieldmap = imgPara.fieldmap;
                disp('Field map input           : True');
            catch
                imgPara2.fieldmap = zeros([size(imgPara2.mask) length(imgPara.fa)]);
                disp('Field map input           : False');
            end
            % check initial phase map
            try
                imgPara2.pini = imgPara.pini;
                disp('Initial phase input       : True');
            catch
                imgPara2.pini = nan(size(imgPara2.mask));
                disp('Initial phase input       : False');
            end
            % check volume fraction of intra-axonal water
            try
                imgPara2.icvf = imgPara.icvf;
                disp('Volume fraction of intra-axonal water input: True');
            catch
                disp('Volume fraction of intra-axonal water input: False');
            end
            % check fibre orientation map
            try
                imgPara2.fo = imgPara.fo;
                disp('Fibre orientation input: True');
            catch
                try
                    imgPara2.theta = imgPara.theta;
                    disp('Fibre orientation input: False');
                    disp('Theta input: True');
                catch
                    if ~algoPara2.DIMWI.isFitFreqMW || ~algoPara2.DIMWI.isFitFreqIW || ~algoPara2.DIMWI.isFitR2sEW
                        error('Fibre orienation map or theta map is required for DIMWI');
                    end
                end
            end
            
            try imgPara2.autoSave = imgPara.autoSave; catch; imgPara2.autoSave = true; end
            disp('Input data is valid.')
            
            % determine the number of estimates
            numEst = 5; % basic setting has 6 estimates
            if obj.epgx_params.isExchange
                numEst = numEst + 1;
            end
            if algoPara2.isFitT1mw
                numEst = numEst + 1;
            end
            if algoPara2.DIMWI.isFitVic
                numEst = numEst + 1;
            end
            if algoPara2.DIMWI.isFitFreqMW
                numEst = numEst + 1;
            end
            if algoPara2.DIMWI.isFitFreqIW
                numEst = numEst + 1;
            end
            if algoPara2.DIMWI.isFitR2sEW
                numEst = numEst + 1;
            end
            if algoPara2.isComplex
                numEst = numEst + size(imgPara2.fieldmap,4) + size(imgPara2.pini,4); % total field and inital phase
            end
            algoPara2.numEst = numEst;
        
        end

        % display fitting algorithm information
        function display_message(obj,algoPara)
        
            %%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
            disp('--------------');
            disp('Fitting option');
            disp('--------------');
            if algoPara.isNormData
                disp('GRE data is normalised before fitting');
            else
                disp('GRE data is not normalised before fitting');
            end
            if algoPara.isFastEPG
                disp('Fast EPG approach is used');
                warning('It does not yield the exact result with default fitting with EPG');
            else
                disp('Fast EPG approach is not used');
            end
            fprintf('Max. iterations    : %i\n',algoPara.maxIter);
            fprintf('Function tolerance : %.2e\n',algoPara.fcnTol);
            fprintf('Step tolerance     : %.2e\n',algoPara.stepTol);
            % type of fitting
            if algoPara.isComplex==0
                disp('Fitting with complex model');
            else 
                disp('Fitting with magnitude model');
            end
            % initial guess and fitting bounds
            if isempty(algoPara.userDefine.x0)
                disp('Default initial guess     : True');
            else
                disp('Default initial guess     : False');
            end
            if isempty(algoPara.userDefine.lb)
                disp('Default lower bound       : True');
            else
                disp('Default lower bound       : False');
            end
            if isempty(algoPara.userDefine.ub)
                disp('Default upper bound       : True');
            else
                disp('Default upper bound       : False');
            end
            % initial guess for in-vivo case
            if algoPara.isInvivo
                disp('Initial guesses for in vivo study');
            else
                disp('Initial guesses for ex vivo study');
            end
            
            disp('Cost function options:');
            if algoPara.isWeighted
                disp('Cost function weighted by echo intensity: True');
                disp(['Weighting method: ' algoPara.weightMethod]);
                if strcmpi(algoPara.weightMethod,'1stEcho')
                    disp(['Weighting power: ' num2str(algoPara.weightPower)]);
                end
            else
                disp('Cost function weighted by echo intensity: False');
            end
            
            disp('--------------------------');
            disp('Multi-compartment T1 model');
            disp('--------------------------');
            if obj.epgx_params.isExchange
                disp('Exchange - to be fitted');
            else
                disp('Exchange - no exchange');
            end
            if algoPara.isFitT1mw
                disp('T1mw - to be fitted');
            else
                disp('T1mw - fixed');
                fprintf('T1mw will be fixed as %.3f s\n',obj.t1_mw);
            end
            if obj.epgx_params.isEPG
                disp('EPG - True');
                fprintf('No. of RF pulses to equilibrium: %i\n', algoPara.npulse);
                fprintf('Initial RF phase               : %i\n', algoPara.rfphase);
            else
                disp('EPG - no EPG simulation');
            end
            
            disp('------------------------------------');
            disp('Diffusion informed MWI model options');
            disp('------------------------------------');
            if ~algoPara.DIMWI.isFitVic
                disp('Intra-axonal and extra-cellular water signal are fitted together');
            else
                disp('Intra-axonal and extra-cellular water signal are fitted separately.');
            end
            if ~algoPara.DIMWI.isFitFreqMW
                disp('Frequency - myelin water estimated by HCFM');
            else
                disp('Frequency - myelin water to be fitted');
            end
            if ~algoPara.DIMWI.isFitFreqIW
                disp('Frequency - intra-axonal water estimated by HCFM');
            else
                disp('Frequency - intra-axonal water to be fitted');
            end
            if ~algoPara.DIMWI.isFitR2sEW
                disp('R2*       - extra-cellular water estimated by HCFM');
            else
                disp('R2*       - extra-cellular water to be fitted');
            end
            
            disp('-------------------------------')
            disp('Parameter to be fixed for DIMWI')
            disp('-------------------------------')
            disp(['Field strength (T)                       : ' num2str(obj.B0)]);
            disp(['B0 direction(x,y,z)                      : ' num2str(obj.B0dir(:)')]);
            disp(['Relative myelin water density            : ' num2str(obj.rho_mw)]);
            disp(['Myelin isotropic susceptibility (ppm)    : ' num2str(obj.x_i)]);
            disp(['Myelin anisotropic susceptibility (ppm)  : ' num2str(obj.x_a)]);
            disp(['Exchange term (ppm)                      : ' num2str(obj.E)]);
            
        end

        % set up DIMWI required data
        function [icvf, ff, theta] = setup_DIMWI(obj,imgParam, algoParam)
            
            DIMWIParams = algoParam.DIMWI;
            
            dims	= size(imgParam.img);   
            dims    = dims(1:3);
            
            % check if intra-axonal water volume fraction is needed
            if ~DIMWIParams.isFitVic
                icvf = double(imgParam.icvf);
            else
                icvf = zeros(dims);
            end
            
            % check if fibre orientation is needed
            if ~DIMWIParams.isFitFreqMW || ~DIMWIParams.isFitFreqIW || ~DIMWIParams.isFitR2sEW
                ff    = double(imgParam.ff); % fibre fraction
                try 
                    fo    = double(imgParam.fo); % fibre orientation w.r.t. B0
                    theta = zeros(size(ff));
                    for kfo = 1:size(fo,5)
                        theta(:,:,:,kfo) = AngleBetweenV1MapAndB0(fo(:,:,:,:,kfo),obj.B0dir);
                    end
                catch 
                    theta = double(imgParam.theta); % theta map
                end
                % normalise fibre fraction
                ff              = bsxfun(@rdivide,ff,sum(ff,4));
                ff(isnan(ff))   = 0;
            else
                ff      = ones(dims);
                theta   = zeros(dims);
            end
            
        end

        % set up initial starting points and fitting boundary
        function [x0,lb,ub] = set_up_x0lbub(obj,y,extraData,DIMWIParams,fitParams)
            % get initial guesses from extraData
            t1iew0  = extraData.t1iew0;
            if isfield(extraData,'mwf0');    mwf0       = min(extraData.mwf0,0.3); end % maximum starting MWF is 30%
            if isfield(extraData,'t2siew0'); r2siew0    = 1/extraData.t2siew0;  end 
            if ~exist("r2siew0",'var');      r2siew0    = extraData.r2s0;       end
            rho0    = max(extraData.m00,max(abs(y(:))));
            db0     = extraData.db0;
            pini0   = extraData.pini0;

            % get user defined bounds
            userDefine = fitParams.userDefine;

            % set initial starting points and fitting boundary
            if fitParams.isInvivo

                % range for in vivo, for 3T
                if fitParams.isSelfBoundary
                    mwi_initial_guess_3T_invivo_selfboundary;
                else
                    mwi_initial_guess_3T_invivo;
                end
            
            else
                % range for ex vivo
                mwi_initial_guess_3T_exvivo
            end

            % volume fraction of intra-axonal water
            if ~DIMWIParams.isFitVic
                Amy0 = mwf0*rho0;            Amylb = 0;	Amyub = 1.5*rho0;
                Aie0 = (1-mwf0)*rho0;        Aielb = 0;	Aieub = 1.5*rho0;
                
                x0 = [Amy0 ,Aie0 ,r2smw0 ,r2siw0 ];
                lb = [Amylb,Aielb,r2smwlb,r2siwlb];
                ub = [Amyub,Aieub,r2smwub,r2siwub];
            else
                Amy0 = mwf0*rho0;            Amylb = 0;	Amyub = 1.5*rho0;
                Aiw0 = iwf0*rho0;            Aiwlb = 0;	Aiwub = 1.5*rho0;
                Aew0 = ewf0*rho0;            Aewlb = 0;	Aewub = 1.5*rho0;
                
                x0 = [Amy0 ,Aiw0 ,Aew0, r2smw0 ,r2siw0 ];
                lb = [Amylb,Aiwlb,Aewlb,r2smwlb,r2siwlb];
                ub = [Amyub,Aiwub,Aewub,r2smwub,r2siwub];
            end
            % R2* of extracellular water
            if DIMWIParams.isFitR2sEW
                x0 = [x0,r2sew0]; lb = [lb,r2sewlb]; ub = [ub,r2sewub];
            end
            % T1 and exchange
            if fitParams.isFitT1mw
                x0 = [x0,t1s0]; lb = [lb,t1slb]; ub = [ub,t1sub];
            end
            x0 = [x0,t1l0]; lb = [lb,t1llb]; ub = [ub,t1lub];
            if obj.epgx_params.isExchange
                x0 = [x0,kls0]; lb = [lb,klslb]; ub = [ub,klsub];
            end
            
            % Frequency of myelin water
            if DIMWIParams.isFitFreqMW
                freq_mw0 = 5 *obj.B0/3;	freq_mwlb = freq_mw0 - 10*obj.B0;	freq_mwub  = freq_mw0 + 10*obj.B0;
                % freq_mw0 = 5*2*pi;	freq_mwlb = -5*obj.B0*2*pi;	freq_mwub  = freq_mw0 + 10*obj.B0*2*pi;
                
                x0 = double([x0,freq_mw0]); lb = double([lb,freq_mwlb]); ub = double([ub,freq_mwub]);
            end
            % Frequency of intra-axonal water
            if DIMWIParams.isFitFreqIW
                % freq_iw0 = -2*2*pi;	freq_iwlb = freq_iw0 - 8*obj.B0*2*pi;	freq_iwub  = freq_iw0 + 8*obj.B0*2*pi;
                freq_iw0 = -2 *obj.B0/3;	freq_iwlb = freq_iw0 - 8*obj.B0;	freq_iwub  = freq_iw0 + 8*obj.B0;
                
                x0 = double([x0,freq_iw0]); lb = double([lb,freq_iwlb]); ub = double([ub,freq_iwub]);
            end
            
            % extra parameters that depended on fitting method
            if fitParams.isComplex
                % freq_bg0 = db0(:).'*2*pi;	freq_bglb = (db0(:).'-8* obj.B0)*2*pi;      freq_bgub  = (db0(:).'+8*obj.B0)*2*pi;
                freq_bg0 = db0(:).';	    freq_bglb = db0(:).'-8* obj.B0;         freq_bgub  = db0(:).'+8*obj.B0;
                pini00 = pini0(:).';        pinilb    = ones(size(pini0))*(-2*pi);  piniub     = ones(size(pini0))*2*pi;
            
                % extra parameters
                x0 = double([x0,freq_bg0, pini00]);
                lb = double([lb,freq_bglb,pinilb]);
                ub = double([ub,freq_bgub,piniub]);
            end
            if ~isempty(userDefine.x0)
                x0(~isnan(userDefine.x0)) = userDefine.x0(~isnan(userDefine.x0));
            end
            if ~isempty(userDefine.lb)
                lb(~isnan(userDefine.lb)) = userDefine.lb(~isnan(userDefine.lb));
            end
            if ~isempty(userDefine.ub)
                ub(~isnan(userDefine.ub)) = userDefine.ub(~isnan(userDefine.ub));
            end

            if fitParams.DEBUG
                x0
            end
        end

    end

    methods(Static)

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
        
        
        % Extra variables
        
        %%% The maximum order varies through the sequence. This can be used to speed up the calculation
        np = length(phi);
        % if not defined, assume want max
        if ~exist('kmax','var')
            kmax = np - 1;
        end
        
        
        %%% Number of states is 6x(kmax +1) -- +1 for the zero order
        N=6*(kmax+1);
        
        
        % Set up matrices for Relaxation and Exchange
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
            
            
            % Main body of gradient echo sequence, loop over TRs
            
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

        % setup output related operations
        function [output_dir,output_filename] = setup_output(imgPara,output_filename)
        %
        % Input
        % --------------
        % imgPara       : structure array contains all image data
        % Output setting
        % ==============
        %   .output_dir      : directory to store final results (default:
        %                      '/current_directory/mwi_results/')
        %   .output_filename : output filename in text string (default:
        %                      'mwi_results.mat')
        %   .identifier      : temporary file identifier, a 8-digit code (optional)
        %
        % Output
        % --------------
        % output_filename : full output filename with extension
        % temp_filename   : full temporary filename
        % identifier      : 8-digit code (in string)
        %
        % Description: setup output related operations
        %
        % Kwok-shing Chan @ DCCN
        % k.chan@donders.ru.nl
        % Date created: 16 Nov 2020
        % Date modified:
        %
        %
            % check if user specified output directory
            if isfield(imgPara,'output_dir')
                output_dir = imgPara.output_dir;
            else
                output_dir    = fullfile(pwd,'mwi_results');
            end
            if ~exist(output_dir,'dir')
                mkdir(output_dir)
            end
            
            % check if user specified output filename
            if isfield(imgPara,'output_filename')
                output_filename = imgPara.output_filename;
                [~,output_filename,~] = fileparts(output_filename);
                output_filename = [output_filename '.mat'];
            end
            output_filename = fullfile(output_dir,output_filename);
            
        end

        function [temp_dir,temp_prefix, identifier] = setup_temp_dir(imgPara,temp_prefix)
        %
        % Input
        % --------------
        % imgPara       : structure array contains all image data
        % Output setting
        % ==============
        %   .output_dir      : directory to store final results (default:
        %                      '/current_directory/mwi_results/')
        %   .output_filename : output filename in text string (default:
        %                      'mwi_results.mat')
        %   .identifier      : temporary file identifier, a 8-digit code (optional)
        %
        % Output
        % --------------
        % output_filename : full output filename with extension
        % temp_filename   : full temporary filename
        % identifier      : 8-digit code (in string)
        %
        % Description: setup output related operations
        %
        % Kwok-shing Chan @ DCCN
        % k.chan@donders.ru.nl
        % Date created: 16 Nov 2020
        % Date modified:
        %
        %
            % check if user specified output directory
            if isfield(imgPara,'output_dir')
                output_dir = imgPara.output_dir;
            else
                output_dir = fullfile(pwd,'mwi_results');
            end
            if ~exist(output_dir,'dir')
                mkdir(output_dir)
            end
            
            % create a new identifier if not provided
            identifier = [];
            % make sure the identifier is unique
            while isempty(identifier) || exist([temp_prefix identifier '.mat'],'file')
                % identifier = num2str(randi(9));
                identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
                % for k = 2:8; identifier = [identifier,num2str(randi(9))]; end
                % identifier = num2str(identifier);
            end
            
            temp_dir        =  fullfile(output_dir,['temporary_' identifier]);
            temp_prefix   = [temp_prefix identifier];
            
            if ~exist(temp_dir,'dir')
                mkdir(temp_dir)
            end
            
            end
    
        
    end
end