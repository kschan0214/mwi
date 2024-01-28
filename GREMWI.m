classdef GREMWI
% This class implements the equations in 
% Chan, K.-S., Marques, J.P., 2020. Multi-compartment relaxometry and diffusion informed myelin water imaging  
% Promises and challenges of new gradient echo myelin water imaging methods. 
% Neuroimage 221, 117159. https://doi.org/10.1016/j.neuroimage.2020.117159
% Chan, K.-S., Chamberland, M., Marques, J.P., 2023. On the performance of multi-compartment relaxometry for 
% myelin water imaging (MCR-MWI)  test-retest repeatability and inter-protocol reproducibility. 
% Neuroimage 266, 119824. https://doi.org/10.1016/j.neuroimage.2022.119824
%
% Two models are included: (1) GRE-MWI and (2) GRE-DIMWI
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 24 January 2024
% Date modified:
%
    properties (Constant)
            gyro = 42.57747892;
    end

    properties (GetAccess = public, SetAccess = protected)
        
        % tissue parameters
        x_i     = -0.1;     % ppm
        x_a     = -0.1;     % ppm
        E       = 0.02;     % ppm
        rho_mw  = 0.42;     % relative water proton density relative to IEW
        % hardware setting
        B0      = 3; % T
        B0dir   = [0;0;1]; % main magnetic field direction with respect to FOV
        te;

    end

    methods
        %% Constructor
        function obj = GREMWI(te,fixed_params)
            obj.te = double(te(:));
            
            % fixed tissue and scanner parameters
            if nargin == 2
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
            end

        end

        %% Data fitting
        % prepare data for MCR-MWI fitting
        function [algoPara,data_obj_all] = data_preparation(obj,algoPara,imgPara)

            disp('=================================================');
            disp('Myelin water imaing: GRE-(DI)MWI Data Preparation');
            disp('=================================================');

            %%%%%%%%%% create directory for temporary results %%%%%%%%%%
            default_output_filename = 'gremwi_results.mat';
            temp_prefix             = 'temp_gremwi_';
            [output_dir, output_filename]       = obj.setup_output(imgPara,default_output_filename);
            [temp_dir,temp_prefix, identifier]  = obj.setup_temp_dir(imgPara,temp_prefix);

            %%%%%%%%%% validate algorithm and image parameters %%%%%%%%%%
            [algoPara,imgPara] = obj.check_and_set_default(algoPara,imgPara);

            %%%%%%%%%% capture all image parameters %%%%%%%%%%
            data  = double(imgPara.img);
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
            [r2s0,~,m00]    = R2star_trapezoidal(data,obj.te);
            mask_valida_r2s = and(~isnan(r2s0),~isinf(r2s0));
            mask_valid_m0   = m00 > 0;
            r2s0(mask_valida_r2s == 0)  = 0;
            r2s0(mask_valid_m0==0)      = 0;
            m00(mask_valida_r2s == 0)  = 0;
            m00(mask_valid_m0==0)      = 0;


            % only works for 3T data
            advancedStarting = algoPara.advancedStarting;
            if strcmpi(advancedStarting,'default') || strcmpi(advancedStarting,'robust')
                fprintf('Estimate starting points using predefined model...')
                
                if obj.B0 > 2.5 && obj.B0 < 3.5 % 3T
                    t2s_pre = [10e-3,60e-3];    % [T2sMW, T2sIEW] in second
                elseif obj.B0 > 1 && obj.B0 < 2 % 1.5T
                    t2s_pre = [10e-3,80e-3];    
                elseif obj.B0 > 6 && obj.B0 < 8 % 7T
                    t2s_pre = [10e-3,40e-3];    
                end
                
                switch advancedStarting
                    case 'default'
                        [m00,mwf0,t2siew0] = superfast_mwi_2m_standard_self(data,obj.te,t2s_pre(1));
                    case 'robust'
                        [m00,mwf0] = superfast_mwi_2m_standard(data,obj.te,t2s_pre);
                end
                m00 = sum(m00,4); % total water
                % also masked out problematic voxels detected by superfast method
                mask_valid_m0 = m00 > 0;
                disp('Completed.')
            end
            
            % final mask for fitting
            mask = and(and(and(mask,mask_valida_r2s),mask_valid_m0),mask_valid_DIMWI);
            
            %%%%%%%%%% Batch processing preparation %%%%%%%%%%
            % set up data_obj for batch processing
            % input data
            data_obj_all = setup_batch_create_data_obj_slice(data,      mask, 'data');
            data_obj_all = setup_batch_create_data_obj_slice(fm0,       mask, 'fm0',    data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(pini0,     mask, 'pini0',  data_obj_all);
            % DIMWI
            data_obj_all = setup_batch_create_data_obj_slice(icvf,      mask, 'icvf',   data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(theta,     mask, 'theta',  data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(ff,        mask, 'ff',     data_obj_all);
            % derived initial guess
            data_obj_all = setup_batch_create_data_obj_slice(m00,       mask, 'm00',	data_obj_all);
            data_obj_all = setup_batch_create_data_obj_slice(r2s0,      mask, 'r2s0',	data_obj_all);
            if exist('mwf0','var');     data_obj_all = setup_batch_create_data_obj_slice(mwf0,      mask, 'mwf0',       data_obj_all); end
            if exist('t2siew0','var');  data_obj_all = setup_batch_create_data_obj_slice(t2siew0,   mask, 't2siew0',	data_obj_all); end
            
            NTotalSamples = numel(mask(mask>0));
            for kbat = 1:numel(data_obj_all)
                data_obj    = data_obj_all(kbat);
                batchNumber = kbat;
                NSamples    = size(data_obj.data,1);
                save(fullfile(temp_dir,strcat(temp_prefix,'_batch',num2str(kbat))),'data_obj','identifier','scaleFactor','batchNumber',...
                    'NTotalSamples','NSamples');
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
            disp('Myelin water imaing: GRE-(DI)MWI model fitting');
            disp('==============================================');

            % display some messages
            obj.display_algorithm_info(algoPara);
    
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
            default_output_filename = 'mwi_gre_results.mat';
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
            jacobianMode    = algoPara.jacobian;
            isAutoSave      = algoPara.autoSave;

            % data fitting method related parameters
            fitParams.isComplex      = algoPara.isComplex;
            fitParams.isWeighted     = algoPara.isWeighted;
            fitParams.weightMethod   = algoPara.weightMethod;
            fitParams.weightPower    = algoPara.weightPower;
            fitParams.userDefine     = algoPara.userDefine;
            fitParams.isInvivo       = algoPara.isInvivo;
            fitParams.isSelfBoundary = algoPara.isSelfBoundary;
            fitParams.DEBUG          = algoPara.DEBUG;

            % DIMWI parameters
            DIMWI = algoPara.DIMWI;

            %%%%%%%%%% set fitting options %%%%%%%%%%
            options = optimoptions(@lsqnonlin,'MaxIter',maxIter,'MaxFunctionEvaluations',200*numEst,...
                'StepTolerance',stepTol,'FunctionTolerance',fcnTol,'Display','off','Jacobian',jacobianMode);

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
                
                %%%%%%%%%% progress display %%%%%%%%%%

                disp('--------')
                disp('Progress')
                disp('--------')

                %%%%%%%%%% fitting main %%%%%%%%%%
                for kbat = 1:Nbatch
                    if isBatch
                        % if batch mode then check if the input mat
                        % contains output variable
                        temp_mat = strcat(temp_files,['_batch' num2str(kbat)]);
                        variableInfo = who('-file', temp_mat);
                        isFit = ~ismember('res_obj', variableInfo);
                    else
                        NSamples    = size(data_obj.data,1);
                        isFit       = true;
                    end

                    % display current batch number
                    fprintf('#Slice %3d/%d',kbat,Nbatch);

                    if isFit 
                        if isBatch
                            % if batch mode then load data from disk
                            load(temp_mat,'data_obj','NSamples');
                        end

                        % display number of voxels to be fitted
                        fprintf(',   #Voxel = %6d',NSamples);

                        % process non-empty slices
                        if ~data_obj.isEmptySlice

                            data    = double(data_obj.data);
                            
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
                                
                                if isfield(data_obj,'mwf0');    extra_data(k).mwf0      = double(data_obj.mwf0(k));         end % MWF
                                if isfield(data_obj,'t2siew0'); extra_data(k).t2siew0   = double(data_obj.t2siew0(k));      end
                            end
                            
                            % start timer
                            start = tic;
                            % create an empty array for fitting results 
                            resnorm   = zeros(size(data,1),1);
                            iter      = zeros(size(data,1),1);
                            exitflag  = zeros(size(data,1),1);
                            if isParallel
                                parfor k = 1:size(data,1)
                                    s           = squeeze(data(k,:,:));
                                    extraData   = extra_data(k);
                                    
                                        [fitRes(k), resnorm(k),exitflag(k),output] = ...
                                            obj.fit(s,extraData,DIMWI,fitParams,options);
                                        iter(k) = output.iterations;
                                end
                            else
                                for k = 1:size(data,1)
                                    s           = squeeze(data(k,:,:));
                                    extraData   = extra_data(k);
                                    
                                        [fitRes(k), resnorm(k),exitflag(k),output] = ...
                                            obj.fit(s,extraData,DIMWI,fitParams,options);
                                        iter(k) = output.iterations;
                                end
                            end
                            
                            % display processing time
                            D = duration(0,0,toc(start),'Format','mm:ss');
                            fprintf(',    Elapsed(MM:SS): %s \n',string(D));

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
                            res_obj.timeElapsed  = D;
                            
                            % display progress
                            if isBatch
                                save(temp_mat,'res_obj','-append')
                            end
                        else
                            fprintf('\n');
                            
                        end
                    else
                        if isBatch
                            % if batch mode then load data from disk
                            load(temp_mat,'NSamples','res_obj');
                            % display number of voxels to be fitted
                            fprintf(',   #Voxel = %6d',NSamples);
                        end
                        fprintf(',    Elapsed(MM:SS): %s ',string(res_obj.timeElapsed ));
                        fprintf(',    this slice is previously processed.\n')
                    end
                end

                if isBatch
                    
                    res_obj = []; isSetup = false;
                    for kbat = 1:Nbatch
                        temp_mat = strcat(temp_files,['_batch' num2str(kbat)]);
                        temp = load(temp_mat);
                        if isfield(temp,'res_obj')
                            if ~isSetup
    
                                fieldname = fieldnames(temp.res_obj);
                                for kfield = 1:numel(fieldname)
                                    if numel(temp.res_obj.(fieldname{kfield})) ~= 1
                                        dims = [size(temp.res_obj.(fieldname{kfield}),1),size(temp.res_obj.(fieldname{kfield}),2),...
                                            Nbatch,size(temp.res_obj.(fieldname{kfield}),3)];
                                        res_obj.(fieldname{kfield}) = zeros(dims);
                                    end
                                end
                                isSetup = true;
                            end
                            fieldname = fieldnames(res_obj);
                            for kfield = 1:numel(fieldname)
                                res_obj.(fieldname{kfield})(:,:,kbat,:) = temp.res_obj.(fieldname{kfield});
                            end
                        end
                    end
                    delete(strcat(temp_files,'_batch*'));
                    rmdir(algoPara.temp_dir);
                end
                if isAutoSave
                    save(output_filename,'res_obj','-append');
                end

                disp('The process is completed.')

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
        function [fitRes, resnorm,exitflag,output] = fit(obj,y,extraData, DIMWI, fitParams,options)

            % if the fitting option is not provided then set up here
            if nargin < 6
                %%%%%%%%%% set fitting options %%%%%%%%%%
                options = optimoptions(@lsqnonlin,'MaxIter',fitParams.maxIter,'MaxFunctionEvaluations',200*fitParams.numEst,...
                    'StepTolerance',fitParams.stepTol,'FunctionTolerance',fitParams.fcnTol,'Display','off',...
                    'Jacobian',fitParams.jacobian);
            end

            %%%%%%%%%% Step 1: Prepare fitting algorithm %%%%%%%%%%
            % set the starting point, lower bound and upper bound
            [x0,lb,ub] = obj.set_up_x0lbub(y,extraData,DIMWI,fitParams);

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
            try
                [x,resnorm,~,exitflag,output] = lsqnonlin(@(x)obj.residuals(x,y,DIMWI, w,fitParams),x0,lb,ub,options);
            catch
                x = zeros(size(x0));
                resnorm = nan;
                exitflag = nan;
                output.iterations = nan;
            end

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
            if DIMWI.isFitFreqMW
                fitRes.Freq_MW = x(counter);    counter= counter + 1;
            end
            if DIMWI.isFitFreqIW
                fitRes.Freq_IW = x(counter);    counter= counter + 1;
            end
            if fitParams.isComplex
                fitRes.Freq_BKG = x(counter);   counter= counter + 1;
                fitRes.pini = x(counter:end);
            end

        end

        % compute the residuals and jacobian for data fitting
        function [residuals, jacobian] = residuals(obj, pars, y, DIMWI, w, fitParams)

            % forward model 
            s = obj.FWD(pars, DIMWI, fitParams);

            if ~fitParams.isComplex
                s = abs(s);
                y = abs(y);
            end

            % compute resiudal (w/ weighting)
            residuals = s(:) - y(:);
            residuals = residuals .* w(:);

            % if complex data fitting then separate real/imaginary
            if fitParams.isComplex
                residuals = [real(residuals);imag(residuals)];
            end

            % Jacobian
            if nargout > 1
                % numerical jacobian
                jacobian = zeros(numel(s),numel(pars));
                for k = 1:numel(pars)
                    h               = 1e-12;
                    pars_h          = pars;
                    pars_h(k)       = pars(k) + h;
                    if fitParams.isComplex
                        tmp             = (obj.FWD(pars_h, DIMWI, fitParams) - s)/h;
                    else
                        tmp             = (abs(obj.FWD(pars_h, DIMWI, fitParams)) - s)/h;
                    end
                    jacobian(:,k)   = tmp(:) .* w(:);
                end
                if fitParams.isComplex
                    jacobian  = [real(jacobian);imag(jacobian)];
                end
            end
                
            
        end

        % MCR-(DI)MWI Forward model fot data fitting
        function s = FWD(obj, pars, DIMWI, fitParams)

            % get parameters
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
            % frequency shifts
            if DIMWI.isFitFreqMW
                freq_mwbg = pars(counter); counter = counter + 1;
            else
                freq_mwbg = 0;
            end
            if DIMWI.isFitFreqIW
                freq_iwbg = pars(counter); counter = counter + 1;
            else
                freq_iwbg = 0;
            end

            % external effects
            if ~fitParams.isComplex % magnitude fitting
                fbg     = 0;                          
                pini    = 0;
            else    % other fittings
                fbg     = pars(counter);  counter = counter + 1;
                pini    = pars(counter);
            end
            
            freq_mw = freq_mwbg - fbg;
            freq_iw = freq_iwbg - fbg;
            freq_ew = 0;

            %%%%%%%%%% simulate signal based on parameter input %%%%%%%%%%
            s = obj.model_GREDIMWI_nFOD(Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI);
            
        end

        %% Signal models
        % GRE-(DI)MWI T2* signal model for 1 fibre direction
        function s = model_GREDIMWI_1FOD(obj,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI)
            
            % Step 1: validate input variables

            if ~DIMWI.isFitFreqMW || ~DIMWI.isFitFreqIW || ~DIMWI.isFitR2sEW
                hcfmObj = HCFM(obj.te,obj.B0);
            end
            
            % if no phase offset set then no phase offset
            if nargin < 12
                pini=0;
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
                if ~DIMWI.isFitFreqIW 
                    freq_iw = hcfmObj.FrequencyAxon(obj.x_a,g,DIMWI.theta);
                end
            end
            
            %%%%%%%%%% extra decay on extracellular water estimated by HCFM %%%%%%%%%%
            if ~DIMWI.isFitR2sEW
                
                % assume extracellular water has the same T2* as intra-axonal water
                t2sew = t2siw;

                fvf = hcfmObj.FibreVolumeFraction(abs(Aiw),abs(Aew),abs(Amw)/obj.rho_mw);
                % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
                d_e = hcfmObj.DephasingExtraaxonal(fvf,g,obj.x_i,obj.x_a,DIMWI.theta);
                d_e = reshape(d_e,size(obj.te));
                
            else
                d_e = zeros(size(obj.te));
            end
            
            %% T2* decay effect
            s = (Amw .* exp(obj.te * (-1/t2smw + 1i*2*pi*freq_mw)) + ...                % S_MW
                 Aiw .* exp(obj.te * (-1/t2siw + 1i*2*pi*freq_iw)) + ...                % S_IW
                 Aew .* exp(obj.te * (-1/t2sew + 1i*2*pi*freq_ew)).*exp(-d_e)) .* ... % S_EW
                 exp(1i*pini) .* exp(1i*2*pi*fbg.*obj.te);                            % phase offset and total field
            
        end

        % GRE-(DI)MWI T2* signal model for all fibre directions
        function s = model_GREDIMWI_nFOD(obj,Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI)
            
            Nte = numel(obj.te);

            DIMWI.ff = DIMWI.ff ./ sum(DIMWI.ff);
            
            s = zeros([Nte length(DIMWI.ff)]);
            DIMWI_curr = DIMWI;
            for kfo = 1:length(DIMWI.theta)

                DIMWI_curr.theta = DIMWI.theta(kfo);
                
                s(:,kfo) = obj.model_GREDIMWI_1FOD(Amw,Aiw,Aew,t2smw,t2siw,t2sew,freq_mw,freq_iw,freq_ew,fbg,pini,DIMWI_curr);
            end
            s = sum(bsxfun(@times,s,permute(DIMWI.ff(:),[2 1])),2);

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
            try algoPara2.jacobian          = algoPara.jacobian;        catch; algoPara2.jacobian = 'off';end
            
            %%%%%%%%%% 2. check data integrity %%%%%%%%%%
            disp('-----------------------');
            disp('Checking data integrity');
            disp('-----------------------');
            % check if the number of echo times matches with the data
            if numel(obj.te) ~= size(imgPara.img,4)
                error('The length of TE does not match with the 4th dimension of the image.');
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
            
            try algoPara2.autoSave = algoPara.autoSave; catch; algoPara2.autoSave = true; end
            disp('Input data is valid.')
            
            % determine the number of estimates
            numEst = 5; % basic setting has 6 estimates
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
                numEst = numEst + 1 + 1; % total field and inital phase
            end
            algoPara2.numEst = numEst;
        
        end

        % display fitting algorithm information
        function display_algorithm_info(obj,algoPara)
        
            %%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
            disp('--------------');
            disp('Fitting option');
            disp('--------------');
            if algoPara.isNormData
                disp('GRE data is normalised before fitting');
            else
                disp('GRE data is not normalised before fitting');
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
            if isfield(extraData,'mwf0');    mwf0       = min(extraData.mwf0,0.4); end % maximum starting MWF is 40%
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
            % Frequency of myelin water
            if DIMWIParams.isFitFreqMW
                freq_mwbg0 = db0;	freq_mwbglb = db0 - 10*obj.B0;	freq_mwbgub  = db0 + 10*obj.B0;
                
                x0 = double([x0,freq_mwbg0]); lb = double([lb,freq_mwbglb]); ub = double([ub,freq_mwbgub]);
            end
            % Frequency of intra-axonal water
            if DIMWIParams.isFitFreqIW
                % freq_iw0 = -2*2*pi;	freq_iwlb = freq_iw0 - 8*obj.B0*2*pi;	freq_iwub  = freq_iw0 + 8*obj.B0*2*pi;
                freq_iwbg0 = db0;	freq_iwbglb = db0 - 8*obj.B0;	freq_iwbgub  = db0 + 8*obj.B0;
                
                x0 = double([x0,freq_iwbg0]); lb = double([lb,freq_iwbglb]); ub = double([ub,freq_iwbgub]);
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

    methods (Static)

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