% run_preproc_mri_stpy

%addpath /home/common/matlab/spm12/
%spm fmri

% Study structure

projectDir = '/project/3017031.05/exp02/';

iostruct = struct('projectDir', projectDir, ...
                  'rawDir', fullfile(projectDir,'data','raw','mri'), ...
                  'analyzedDir', fullfile(projectDir,'data','analyzed'), ...
                  'codeDir', fullfile(projectDir,'src','code','analysis','preprocessing'), ...
                  'batchDir', fullfile(projectDir,'src','code','analysis','preprocessing','batches'), ...
                  'searchPaths', {''});
              

methods = struct('MECombination',struct('method','paid', ...
                                        'TE',[8.5, 19.3, 30, 41]));
              
% Fill in subject data
% -------------------------------------------------------------------------

% % First subject is subj00. Matlab starts counting at 1, so do subj00
% % separately. In preproc_mri_stpy.m, for this one occasion, search for
% % "iSubject" and change the occurences that are NOT in a call to
% % settings.subj(iSubject) to "(iSubject - 1)" (except for the first occurance
% % where iSubject is defined).
% % .........................................................................
% subj(1).data.dicomPrefix    = 'BRAZAN_20160523_DS00.MR.BRAZAN_SKYRA';
% subj(1).data.anat.ix        = {[15]};
% subj(1).data.anat.name      = {'t1'};
% subj(1).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
% subj(1).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};
% % .........................................................................

subj(1).data.dicomPrefix    = 'BRAZAN_20160525_SUBJ01.MR.BRAZAN_SKYRA';
subj(1).data.anat.ix        = {[15]};
subj(1).data.anat.name      = {'t1'};
subj(1).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(1).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(2).data.dicomPrefix    = 'BRAZAN_20160526_SUBJ02.MR.BRAZAN_SKYRA';
subj(2).data.anat.ix        = {[15]};
subj(2).data.anat.name      = {'t1'};
subj(2).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(2).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(3).data.dicomPrefix    = 'BRAZAN_20160601_SUBJ03.MR.BRAZAN_SKYRA';
subj(3).data.anat.ix        = {[15]};
subj(3).data.anat.name      = {'t1'};
subj(3).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(3).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(5).data.dicomPrefix    = 'BRAZAN_20160526_SUBJ05.MR.BRAZAN_SKYRA';
subj(5).data.anat.ix        = {[15]};
subj(5).data.anat.name      = {'t1'};
subj(5).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(5).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(6).data.dicomPrefix    = 'BRAZAN_20160602_SUBJ06.MR.BRAZAN_SKYRA';
subj(6).data.anat.ix        = {[15]};
subj(6).data.anat.name      = {'t1'};
subj(6).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(6).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(7).data.dicomPrefix    = 'BRAZAN_20160606_SUBJ07.MR.BRAZAN_SKYRA';
subj(7).data.anat.ix        = {[15]};
subj(7).data.anat.name      = {'t1'};
subj(7).data.func.ix        = {[11:14],[16:19],[7:10]};
subj(7).data.func.name      = {'stpy_run1','stpy_run2','prescan'};

subj(8).data.dicomPrefix    = 'BRAZAN_20160530_SUBJ08.MR.BRAZAN_SKYRA';
subj(8).data.anat.ix        = {[15]};
subj(8).data.anat.name      = {'t1'};
subj(8).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(8).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(9).data.dicomPrefix    = 'BRAZAN_20160606_SUBJ09.MR.BRAZAN_SKYRA';
subj(9).data.anat.ix        = {[15]};
subj(9).data.anat.name      = {'t1'};
subj(9).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(9).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(10).data.dicomPrefix    = 'BRAZAN_20160607_SUBJ10.MR.BRAZAN_SKYRA';
subj(10).data.anat.ix        = {[15]};
subj(10).data.anat.name      = {'t1'};
subj(10).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(10).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(11).data.dicomPrefix    = 'BRAZAN_20160606_SUBJ11.MR.BRAZAN_SKYRA';
subj(11).data.anat.ix        = {[15]};
subj(11).data.anat.name      = {'t1'};
subj(11).data.func.ix        = {[11:14],[16:19],[7:10]};
subj(11).data.func.name      = {'stpy_run1','stpy_run2','prescan'};

subj(12).data.dicomPrefix    = 'BRAZAN_20160609_SUBJ12.MR.BRAZAN_SKYRA';
subj(12).data.anat.ix        = {[15]};
subj(12).data.anat.name      = {'t1'};
subj(12).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(12).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(13).data.dicomPrefix    = 'BRAZAN_20160608_SUBJ13.MR.BRAZAN_SKYRA';
subj(13).data.anat.ix        = {[15]};
subj(13).data.anat.name      = {'t1'};
subj(13).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(13).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(14).data.dicomPrefix    = 'BRAZAN_20160613_SUBJ14.MR.BRAZAN_SKYRA';
subj(14).data.anat.ix        = {[15]};
subj(14).data.anat.name      = {'t1'};
subj(14).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(14).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(15).data.dicomPrefix    = 'BRAZAN_20160613_SUBJ15.MR.BRAZAN_SKYRA';
subj(15).data.anat.ix        = {[15]};
subj(15).data.anat.name      = {'t1'};
subj(15).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(15).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(16).data.dicomPrefix    = 'BRAZAN_20160620_SUBJ16.MR.BRAZAN_SKYRA';
subj(16).data.anat.ix        = {[15]};
subj(16).data.anat.name      = {'t1'};
subj(16).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(16).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(17).data.dicomPrefix    = 'BRAZAN_20160621_SUBJ17.MR.BRAZAN_SKYRA';
subj(17).data.anat.ix        = {[15]};
subj(17).data.anat.name      = {'t1'};
subj(17).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(17).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(18).data.dicomPrefix    = 'BRAZAN_20160620_SUBJ18.MR.BRAZAN_SKYRA';
subj(18).data.anat.ix        = {[15]};
subj(18).data.anat.name      = {'t1'};
subj(18).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(18).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(19).data.dicomPrefix    = 'BRAZAN_SUBJ19_20160621.MR.BRAZAN_SKYRA';
subj(19).data.anat.ix        = {[15]};
subj(19).data.anat.name      = {'t1'};
subj(19).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(19).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(20).data.dicomPrefix    = 'BRAZAN_20160629_SUBJ20.MR.BRAZAN_SKYRA';
subj(20).data.anat.ix        = {[15]};
subj(20).data.anat.name      = {'t1'};
subj(20).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(20).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(21).data.dicomPrefix    = 'BRAZAN_20160629_SUBJ21.MR.BRAZAN_SKYRA';
subj(21).data.anat.ix        = {[15]};
subj(21).data.anat.name      = {'t1'};
subj(21).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(21).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(22).data.dicomPrefix    = 'BRAZAN_20160705_SUBJ22.MR.BRAZAN_SKYRA';
subj(22).data.anat.ix        = {[15]};
subj(22).data.anat.name      = {'t1'};
subj(22).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(22).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(23).data.dicomPrefix    = 'BRAZAN_20160705_SUBJ23.MR.BRAZAN_SKYRA';
subj(23).data.anat.ix        = {[15]};
subj(23).data.anat.name      = {'t1'};
subj(23).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(23).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

subj(24).data.dicomPrefix    = 'BRAZAN_20160705_SUBJ24.MR.BRAZAN_SKYRA';
subj(24).data.anat.ix        = {[15]};
subj(24).data.anat.name      = {'t1'};
subj(24).data.func.ix        = {[11:14],[16:19],[7:10],[20:23]};
subj(24).data.func.name      = {'stpy_run1','stpy_run2','prescan','rest'};

% Combine everything in settings structure
settings = struct('io',iostruct, ...
                  'methods',methods, ...
                  'subj',subj);

subjectIx = [1];

% Start preprocessing
preproc_mri_stpy(settings,subjectIx)