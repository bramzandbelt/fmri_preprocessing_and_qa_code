function preproc_mri_stpy(settings,subjectIx)
%PREPROC_MRI_STPY Imaging pre-processing pipeline of selective stopping study
%
% DESCRIPTION
% Preprocessing includes the following steps
% - DICOM to NIfTI conversion
% - MultiEcho combination
% - 
%
%
% INPUTS
% - settings        :
% - subjectIx       : vector of integers, corresponding to subject indices
% 
% DEPENDENCIES
%
%
% Bram Zandbelt, June 2016
% 

spm('defaults', 'FMRI');


todo.dicomConvert           = true;
todo.brainExtraction        = true;
todo.multiEchoCombine       = true;
todo.coregNormalizeSmooth   = true;
todo.removeObsoleteFiles    = false;


for ix = 1:numel(subjectIx)
    
    iSubject = subjectIx(ix);
    
    % Subject directories
    % ---------------------------------------------------------------------
    
    dirs.raw        = fullfile(settings.io.rawDir, sprintf('subj%.2d',iSubject));
    dirs.analyzed   = fullfile(settings.io.analyzedDir, sprintf('subj%.2d',iSubject));
    dirs.jobs       = fullfile(settings.io.analyzedDir, sprintf('subj%.2d',iSubject),'jobs');
    dirs.reports    = fullfile(settings.io.analyzedDir, sprintf('subj%.2d',iSubject),'reports');
    dirs.anat       = fullfile(settings.io.analyzedDir, sprintf('subj%.2d',iSubject),'anat');
    dirs.func       = fullfile(settings.io.analyzedDir, sprintf('subj%.2d',iSubject),'func');
    
    % Subject data
    % ---------------------------------------------------------------------
    anatIx          = settings.subj(iSubject).data.anat.ix;
    anatName        = settings.subj(iSubject).data.anat.name;
    funcIx          = settings.subj(iSubject).data.func.ix;
    funcName        = settings.subj(iSubject).data.func.name;
    
    % Make dirs if they don't exist
    dirNames = fieldnames(dirs);
    for iDir = 1:numel(dirNames)
        if ~exist(getfield(dirs,dirNames{iDir}),'dir')
            mkdir(getfield(dirs,dirNames{iDir}));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. CONVERT DICOM TO NIFTI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if todo.dicomConvert
    
        % Load job
        % =================================================================
        moduleIx    = 1;
        nModule     = numel(settings.subj(iSubject).data.anat.ix) +  ...
                      numel(settings.subj(iSubject).data.func.ix);

        jobTemplate = fullfile(settings.io.batchDir,'dicom_import.m');
        jobId       = cfg_util('initjob', jobTemplate);
        jobs        = repmat(jobTemplate, 1, nModule);
        inputs      = cell(2, nModule);

        modJobId    = 1;

        for i = 1:(nModule - 1)
            cfg_util('replicate', jobId,i);
        end

        % Define inputs
        % =================================================================

        % Anatomical scans
        % -----------------------------------------------------------------
        
        for iRun = 1:numel(anatIx)

            % Select files
            % .............................................................
            anatFiles = '';
            filt = sprintf('^%s.%.4d.*IMA$', ...
                           settings.subj(iSubject).data.dicomPrefix, ...
                           anatIx{iRun});
            anatFiles   = strvcat(anatFiles, ...
                                  spm_select('FPList',dirs.raw,filt));

            % Define output directory
            % .............................................................
            outDir = fullfile(dirs.analyzed,'anat',anatName{iRun});

            if ~ exist(outDir,'dir')
                mkdir(outDir);
            end

            % Fill in job
            % .............................................................
            % 1. DICOM files
            % 2. Output directory

            inputs{1, moduleIx} = cellstr(spm_select('FPList',dirs.raw,filt));
            inputs{2, moduleIx} = cellstr(outDir);

            % Move to next module
            % .............................................................
            moduleIx = moduleIx + 1;

        end

        % Functional scans
        % -----------------------------------------------------------------
        
        for iRun = 1:numel(settings.subj(iSubject).data.func.ix)

            % Select files
            % .............................................................
            funcFiles = '';

            for iEcho = 1:numel(funcIx{iRun})
                filt        = sprintf('^%s.%.4d.*IMA$', ...
                                      settings.subj(iSubject).data.dicomPrefix, ...
                                      funcIx{iRun}(iEcho));
                funcFiles   = strvcat(funcFiles, ...
                                      spm_select('FPList',dirs.raw,filt));

            end

            % Define output directory
            % .............................................................
            outDir = fullfile(dirs.analyzed,'func',funcName{iRun});

            if ~ exist(outDir,'dir')
                mkdir(outDir);
            end

            % Fill in batch
            % .............................................................
            % 1. DICOM files
            % 2. Output directory

            inputs{1, moduleIx} = cellstr(funcFiles);
            inputs{2, moduleIx} = cellstr(outDir);

            % Move to next module
            % .............................................................
            moduleIx = moduleIx + 1;

        end

        % Save and run job
        % =================================================================
        [~,nm,ext] = fileparts(jobTemplate);
        jobFile = fullfile(dirs.jobs,[nm,'_',datestr(now,'yyyymmddTHHMMSS'),ext]);

        sts = cfg_util('filljob', jobId, inputs{:});
        if sts
            cfg_util('savejob', jobId, jobFile);
            
            % Go to report dir, so any graphical output will be written there
            cd(dirs.reports)
            
            spm_jobman('run', jobFile, inputs{:});
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. BRAIN EXTRACTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if todo.brainExtraction
        
        for iRun = 1:numel(anatIx)

            % Select file
            % .............................................................
            anatDir = fullfile(dirs.analyzed,'anat',anatName{iRun});
            filt = sprintf('^s.*nii$');
            anatFile = spm_select('FPList',anatDir,filt);
            
            % Define output file
            % .............................................................
            betFile = fullfile(anatDir,sprintf('anatomy_s%.2d_brain',iSubject));
            
            % Run brain extraction with the FSL Brain Extraction Tool (BET)
            % .............................................................
            system(['bet ' anatFile ' ' betFile '-R -m'])
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. COMBINE MULTI-ECHO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if todo.multiEchoCombine
    
        % 3.1. Realign : Estimate (TE1)
        % =================================================================

        % Load job
        % ----------------------------------------------------------------- 
        moduleIx    = 1;
        nRun        = numel(settings.subj(iSubject).data.func.ix);

        jobTemplate = fullfile(settings.io.batchDir,'realign_estimate.m');
        jobId       = cfg_util('initjob', jobTemplate);
        
        for i = 1:(nRun- 1)
            cfg_util('replicate', jobId,1,i,-1);
        end

        % Select first TE scans of all functional runs
        % ----------------------------------------------------------------- 

        % Replicate session
        for i = 1:(nRun - 1)
            item_mod_id(1) = struct('type','.','subs','val');
            item_mod_id(2) = struct('type','{}','subs',[]);
            item_mod_id(2).subs = num2cell(1);
            cfg_util('setval',jobId,1,item_mod_id,[1 i+1]);
        end
        
        filt        = '.*01.nii$';
        
        for iRun = 1:numel(settings.subj(iSubject).data.func.ix)

            % Select files
            % .............................................................
            funcFiles   = spm_select('FPList',fullfile(dirs.func,funcName{iRun}),filt);

            % Fill in batch
            % .............................................................
            % i. Realign: Estimate: Session i
            
            inputs{iRun, moduleIx} = cellstr(funcFiles);
            
        end

        % Save and run realignment job
        % =================================================================
        [~,nm,ext] = fileparts(jobTemplate);
        jobFile = fullfile(dirs.jobs,[nm,'_',datestr(now,'yyyymmddTHHMMSS'),ext]);

        sts = cfg_util('filljob', jobId, inputs{:});
        if sts
            cfg_util('savejob', jobId, jobFile);
            
            % Go to report dir, so any graphical output will be written there
            cd(dirs.reports)
            
            spm_jobman('run', jobFile, inputs{:});
        end
    
        % 3.2. Apply transformations from TE1 to other TE volumes
        % =================================================================
    
        fprintf(1,[repmat('-',1,72),'\n'])
        fprintf(1,'Running ''Apply realignment parameters of TE1 to other TE volumes'' \n')
        
        for iRun = 1:numel(settings.subj(iSubject).data.func.ix)
            
            funcFiles   = cell(1,numel(funcIx{iRun}));
            
            % Select volumes and group according to TE
            % -------------------------------------------------------------
            for iEcho = 1:numel(funcIx{iRun})
                filt = sprintf('.*%.2d.nii$',iEcho);
                funcFiles{iEcho} = spm_select('FPList',fullfile(dirs.func,funcName{iRun}),filt);
            end
            
            funcFiles = cellfun(@(in1) cellstr(in1), funcFiles, 'Uni',0);
            
            % Assert equal number of volumes for each TE
            assert(all(cell2mat(cellfun(@(in1) size(in1,1),funcFiles,'Uni',0))));
            
            % Apply transformation matrix of TE1 to other TE volumes
            % -------------------------------------------------------------
            nFile = size(funcFiles{1},1);
            
            spm_progress_bar('Init',nFile, ... 
                             sprintf('Applying realignment parameters of TE1 to other TE volumes for %s ...',settings.subj(iSubject).data.func.name{iRun}), ...
                             'Volumes Complete');
            for iFile = 1:nFile
                
                % Read header of TE1 volume
                V       = cell(1,numel(funcIx{iRun}));
                V{1}    = spm_vol(funcFiles{1}{iFile});
                
                % For remaining TE volumes
                for iEcho = 2:numel(funcIx{iRun})
                    
                    % Read header and the volume
                    V{iEcho}        = spm_vol(funcFiles{iEcho}{iFile});
                    I               = spm_read_vols(V{iEcho});
                    
                    % Set transformation matrix equal to TE1 volume
                    V{iEcho}.mat    = V{1}.mat;
                    
                    % Write changes to disk
                    spm_write_vol(V{iEcho},I);
                end
                
                spm_progress_bar('Set',iFile);
            end
            spm_progress_bar('Clear');
        end
        
        % 3.3. Realign : Reslice
        % =================================================================
        %
        %
        %
        
        inputs      = cell(1,1);
        
        moduleIx    = 1;
        nRun        = numel(settings.subj(iSubject).data.func.ix);

        jobTemplate = fullfile(settings.io.batchDir,'realign_reslice.m');
        jobId       = cfg_util('initjob', jobTemplate);
        
        % Select files
        % -----------------------------------------------------------------
        filt        = '.*.nii$';
        funcFiles   = '';
        
        for iRun = 1:numel(settings.subj(iSubject).data.func.ix)

            % Select files
            % .............................................................
            funcFiles   = strvcat(funcFiles, ...
                                  spm_select('FPList',fullfile(dirs.func,funcName{iRun}),filt));

        end
        
        % Fill in batch
        % -----------------------------------------------------------------
        % 1. Realign: Reslice: Images 
        inputs{1, moduleIx} = cellstr(funcFiles);
        
        
        % Save and run Realign: Reslice
        % -----------------------------------------------------------------
        [~,nm,ext] = fileparts(jobTemplate);
        jobFile = fullfile(dirs.jobs,[nm,'_',datestr(now,'yyyymmddTHHMMSS'),ext]);

        sts = cfg_util('filljob', jobId, inputs{:});
        if sts
            cfg_util('savejob', jobId, jobFile);
            
            % Go to report dir, so any graphical output will be written there
            cd(dirs.reports)
            
            spm_jobman('run', jobFile, inputs{:});
        end
        
        % 3.4. Combine multiple echo volumes into one volume
        % =================================================================
        
        MECombineMethod = settings.methods.MECombination.method;
        TE              = settings.methods.MECombination.TE;
        
        
        % Identify which run contains the prescans
        iPrescanRun = find(strncmp(settings.subj(iSubject).data.func.name,'prescan',length('prescan')));
        
        % Select realigned and resliced prescans
        funcFiles   = cell(1,numel(funcIx{iPrescanRun}));
            
        % Step 1 - Select prescan volumes and group according to TE
        % -----------------------------------------------------------------
        for iEcho = 1:numel(funcIx{iPrescanRun})
            filt = sprintf('^r.*%.2d.nii$',iEcho);
            funcFiles{iEcho} = spm_select('FPList',fullfile(dirs.func,funcName{iPrescanRun}),filt);
        end
        
        funcFiles = cellfun(@(in1) cellstr(in1), funcFiles, 'Uni',0);
        
        nFile = size(funcFiles{1},1);
        
        % Step 2 - Read prescan headers and image data
        % -----------------------------------------------------------------
        
        % Get image dimensions of one of the images for pre-allocation
        % purposes
        V1_1 = spm_vol(funcFiles{1}{1});
        
        % Pre-allocate matrices for headers and image data
        V   = cell(1,numel(funcIx{iPrescanRun}));
        Y   = nan([V1_1.dim,nFile,numel(funcIx{iPrescanRun})]);
        
        % Read pre-scan headers and image data
        for iFile = 1:nFile
            for iEcho = 1:numel(funcIx{iPrescanRun})
                V{iEcho} = spm_vol(funcFiles{iEcho}{iFile});
                Y(:,:,:,iFile,iEcho) = spm_read_vols(V{iEcho});
            end
        end
        
        % Step 3 - Compute the weighting of each TE
        % -----------------------------------------------------------------
        switch lower(MECombineMethod)
            case 'paid'
                
                % For each TE, compute temporal signal-to-noise ratio
                % and contrast-to-noise ratio
                for iEcho=1:numel(funcIx{iPrescanRun})
                    tSNR(:,:,:,iEcho) = mean(Y(:,:,:,:,iEcho),4)./std(Y(:,:,:,:,iEcho),0,4);
                    CNR(:,:,:,iEcho) = tSNR(:,:,:,iEcho) * TE(1,iEcho); %% assuming all runs have the same TEs!!
                end

                % Sum contrast to noise ratios across TEs
                CNRTotal = sum(CNR,4);

                % Determine the weighting
                for iEcho=1:size(TE,2)
                    weight(:,:,:,iEcho) = CNR(:,:,:,iEcho) ./ CNRTotal;
                end
                
            case 'te'
                for iEcho=1:size(TE,2)
                    weight(:,:,:,i) = TE(i)/sum(TE);
                end
        end
        
        % Step 4 - Write the weighting images (for reference only)
        % -----------------------------------------------------------------
        
        for iEcho=1:size(TE,2)
            
            VWeights{iEcho} = V{iEcho};
            VWeights{iEcho}.dt = [spm_type('float64'),0];
            
            % New file name
            % .............................................................
            [pth,nme,ext] = fileparts(VWeights{iEcho}.fname);
            VWeights{iEcho}.fname = fullfile(pth,[sprintf('MultiEchoWeights_%smethod_TE%.2d',upper(MECombineMethod),iEcho),ext]);
            
            % Write header of the multi-echo combined volume
            % .........................................................
            spm_create_vol(VWeights{iEcho});
                
            % Write image data of the multi-echo combined volume
            % .........................................................
            spm_write_vol(VWeights{iEcho},weight(:,:,:,iEcho));
            
        end
        
        % Step 5 - Apply the weights to the other functional runs
        % -----------------------------------------------------------------
        
        allFuncRunIx = 1:numel(settings.subj(iSubject).data.func.ix);
        
        for iRun = setdiff(allFuncRunIx,iPrescanRun)
            
            % Select other functional volumes and group according to TE
            % .............................................................
            funcFiles   = cell(1,numel(funcIx{iRun}));
            
            for iEcho = 1:numel(funcIx{iRun})
                filt = sprintf('^r.*%.2d.nii$',iEcho);
                funcFiles{iEcho} = spm_select('FPList',fullfile(dirs.func,funcName{iRun}),filt);
            end
            
            funcFiles = cellfun(@(in1) cellstr(in1), funcFiles, 'Uni',0);
            
            % Assert equal number of volumes for each TE
            assert(all(cell2mat(cellfun(@(in1) size(in1,1),funcFiles,'Uni',0))));
            
            nFile = size(funcFiles{1},1);

            for iFile = 1:nFile
                
                V = cell(1,numel(funcIx{iRun}));
                
                % Read functional volumes headers for all echos
                % .........................................................
                for iEcho = 1:numel(funcIx{iRun})
                    V{iEcho} = spm_vol(funcFiles{iEcho}{iFile});
                end
                
                % Make new filename for combined functional volumes
                % .........................................................
                [pth, nme, ext] = fileparts(funcFiles{1}{iFile});
                
                % combine (c) prefix is added, echo index suffix is removed
                newName         = ['c',nme(1:(max(strfind(nme,'-')-1)))];
                
                % Create header of the multi-echo combined volume
                % .........................................................
                VNew = V{1};
                VNew.fname = fullfile(pth,[newName,ext]);
                
                % Do the actual weighting
                % .........................................................
                YWeighted = zeros(VNew.dim);
                
                for iEcho = 1:numel(funcIx{iRun})
                    Y(:,:,:,iEcho) = spm_read_vols(V{iEcho});
                    YWeighted = YWeighted + Y(:,:,:,iEcho) .* weight(:,:,:,iEcho); 
                end        
                
                % Write header of the multi-echo combined volume
                % .........................................................
                spm_create_vol(VNew);
                
                % Write image data of the multi-echo combined volume
                % .........................................................
                spm_write_vol(VNew,YWeighted);
                
            end
        end        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COREGISTER, NORMALIZATION, AND SMOOTHING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if todo.coregNormalizeSmooth
        
        % Load job
        % ----------------------------------------------------------------- 
        jobTemplate = fullfile(settings.io.batchDir,'coreg_normalize_smooth.m');
        jobId       = cfg_util('initjob', jobTemplate);
        inputs      = cell(3, 1);
        
        
        % Select files
        % ----------------------------------------------------------------- 
        
        % Reference image is the mean functional image, which is located in the directory of first functional run
        filt            = '^mean.*nii$';
        meanFuncImg     = spm_select('FPList',fullfile(dirs.func,funcName{1}),filt);
        
        % Source image is the anatomical scan
        iT1             = find(strncmp(settings.subj(iSubject).data.anat.name,'t1',length('t1')));
        filt            = '.*nii$';
        t1Img           = spm_select('FPList',fullfile(dirs.anat,anatName{iT1}),filt);
        
        % Normalize: Write images are the combined and realigned functional
        % scans
        % Identify which run contains the prescans
        iPrescanRun     = find(strncmp(settings.subj(iSubject).data.func.name,'prescan',length('prescan')));
        allFuncRunIx    = 1:numel(settings.subj(iSubject).data.func.ix);
        funcFiles       = '';
        
        for iRun = setdiff(allFuncRunIx,iPrescanRun)
            filt        = '^cr.*.nii$';
            funcFiles   = strvcat(funcFiles, ...
                                  spm_select('FPList',fullfile(dirs.func,funcName{iRun}),filt));
        end
        
        % Fill in job
        % ----------------------------------------------------------------- 
        % 1. Coregister: Estimate: Reference Image - cfg_files
        % 2. Coregister: Estimate: Source Image - cfg_files
        % 3. Normalise: Write: Images to Write - cfg_files
        
        inputs{1, 1} = cellstr(meanFuncImg);
        inputs{2, 1} = cellstr(t1Img);
        inputs{3, 1} = cellstr(funcFiles);
        
        % -----------------------------------------------------------------
        [~,nm,ext] = fileparts(jobTemplate);
        jobFile = fullfile(dirs.jobs,[nm,'_',datestr(now,'yyyymmddTHHMMSS'),ext]);

        sts = cfg_util('filljob', jobId, inputs{:});
        if sts
            cfg_util('savejob', jobId, jobFile);
            
            % Go to report dir, so any graphical output will be written there
            cd(dirs.reports)
            
            spm_jobman('run', jobFile, inputs{:});
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REMOVE OBSOLETE FILES CREATED IN PREPROCESSING STEPS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if todo.removeObsoleteFiles
        
        fprintf('Removing obsolete preprocessing files of subject %.2d\n',iSubject)
        
        % Select the files to be removed from eacht MRI run (i.e. each
        % directory)
        % =================================================================
        % Non-realigned prescan files to be removed
        % -----------------------------------------------------------------
        filt = '^f.*nii$';
        rmPrescan = spm_select('FPList',fullfile(dirs.func,'prescan'),filt);
        
        % Functional scans: leave only the realigned, combined, normalized
        % and smoothed images (prefix = swcrf) and the txt file containing
        % the realignment parameters
        
        % Resting state scan files to be removed
        % -----------------------------------------------------------------
        filt = '^[^s].*nii$';
        rmRest = spm_select('FPList',fullfile(dirs.func,'rest'),filt);
        
        % fMRI run 1 files to be removed
        % -----------------------------------------------------------------
        filt = '^[^s].*nii$';
        rmRun1 = spm_select('FPList',fullfile(dirs.func,'stpy_run1'),filt);
        
        % fMRI run 2 files to be removed
        % -----------------------------------------------------------------
        filt = '^[^s].*nii$';
        rmRun2 = spm_select('FPList',fullfile(dirs.func,'stpy_run2'),filt);
        
        % The file list (full paths) with all files to be removed
        % =================================================================
        toRemove = cellstr(strvcat(rmPrescan,rmRest,rmRun1,rmRun2));
        
        % Delete all files selected for removal
        % =================================================================
        delete(toRemove{:})
        
    end
end
    
    
    
end
