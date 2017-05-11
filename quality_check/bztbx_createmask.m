function bztbx_createmask(varargin)
% ________________________________________________________________________
% CREATE MASK Construction of binary images based on  
%
% FORMAT: bztbx_createmask('sphere',coord,radius,spacedefimg,tgtdir,outputtype)
% DESCRIPTION       Makes binary image(s) of one or more spheres
% coord             - list of world-space coordinates (N x 3 vector)
% radius            - radius of sphere, in mm (scalar)
% spacedefimg       - filename of space defining image (char array)
% tgtdir            - filename of target directory (char array)
% outputtype        - 'separate','cluster' (number labeled),'mask' (all ones)
%
% FORMAT: bztbx_createmask('box',coord,width,spacedefimg,tgtdir,outputtype)
% DESCRIPTION       Makes binary image(s) of one or more boxes
% dimensions        - dimensions, in mm (1 x 3 vector)
%
% FORMAT: bztbx_createmask('spm',statdir,statspec,spacedefimg,tgtdir,outputtype)
% DESCRIPTION       Makes binary image(s) of activated clusters in SPM
% statdir           - statistics directory (char array)
% statspec          - cell array containing info re:
%                       - contrast number
%                       - significance level (e.g. .05)
%                       - multiple comparison type ('none','FDR','FWE')
%                       - cluster extent threshold
%                       eg. statspec = {1,.05,'FWE',0};  
%
% FORMAT: bztbx_createmask('image',srcfls,spacedefimg,tgtdir,outputtype)
% DESCRIPTION        
% srcfls            - source images
%
% FORMAT: bztbx_createmask('LPBA_mask',srcfls,tgtdir,k)
% DESCRIPTION        
% srcfls            - source images, these should be masks
% k                 - minimum cluster size in new map (optional)
%
% ________________________________________________________________________
% Bram B. Zandbelt, October 2009
% b.b.zandbelt@umcutrecht.nl
%
% version 1.0

% Bugs:
% When making ROIs with outputtype 'separate' matrices contain values like
% 1.000 instead of 1. This is not the case when using cluster or mask.
%

if nargin == 0
    error('Wrong number of input arguments. Type ''help bztbx_createmask'' for further assistance.');
end


switch lower(varargin{1})
    case {'sphere','box'}
        coord = varargin{2};
        dims = varargin{3};
        spacedefimg = varargin{4};
        tgtdir = varargin{5};
        outputtype = varargin{6};
    case 'spm'
        addpath /mnt/home/bram/m-files/thomas_nichols/      % For CorrClusTh function
        statdir = varargin{2};
        if isempty(varargin{3})
             load(fullfile(statdir,'SPM.mat'));
             k = CorrClusTh(SPM,0.001,0.05);
             spc = {1,0.001,'none',k}
        else
            spc = varargin{3};
        end
        spacedefimg = varargin{4};
        tgtdir = varargin{5};
        outputtype = varargin{6};
    case 'image'
        srcfls = varargin{2};
        spacedefimg = varargin{3};
        tgtdir = varargin{4};
        outputtype = varargin{5};
    case 'lpba_mask'
        srcfls = varargin{2};
        tgtdir = varargin{3};
        if nargin == 4
            k = varargin{4};
        else
            k = 0;
        end
end
        

addpath /raid/bram/SPMTOOLBOX/spm5/toolbox/marsbar/                             % Add Marsbar path
marsbar('on');                                                     

% Create VOIs
switch lower(varargin{1})
    case 'sphere'
        
        for ivoi = 1:size(coord,1)
            c = coord(ivoi,:);                                                  % VOI center coordinate
            r = dims;                                                           % VOI radius
            d = sprintf('%0.1fmm radius sphere at [%0.1f %0.1f %0.1f]',r,c);    % VOI description
            l = sprintf('voi_%02.0f_sphere_%0.0f-%0.0f_%0.0f_%0.0f',ivoi,r,c);  % VOI label                % VOI label
            o = maroi_sphere(struct('centre',c,'radius',r));                    % Create VOI object
            
            o = descrip(o,d);                                                           % Add VOI description to VOI object
            o = label(o,l);                                                             % Add VOI label to VOI object
 
            fn = fullfile(tgtdir,[label(o),'_roi.mat']);
            saveroi(o,fn);
            
        end
    case 'box'
        
        for ivoi = 1:size(coord,1)
            c = coord(ivoi,:);                                                          % VOI center coordinate
            w = dims;                                                                   % VOI dimensions
            d = sprintf('[%0.1f %0.1f %0.1f] box at [%0.1f %0.1f %0.1f]',w,c);          % VOI description
            l = sprintf('voi_%02.0f_box_w-%0.0f_%0.0f_%0.0f-%0.0f_%0.0f_%0.0f',ivoi,w,c); % VOI label      % VOI label
            o = maroi_box(struct('centre',c,'widths',w));                               % Create VOI object
            
            o = descrip(o,d);                                                           % Add VOI description to VOI object
            o = label(o,l);                                                             % Add VOI label to VOI object
 
            fn = fullfile(tgtdir,[label(o),'_roi.mat']);
            saveroi(o,fn);
            
        end
        
    case 'spm'
        
        % Load the RFX SPM.mat
        load(fullfile(statdir,'SPM.mat'));
        
        % Create the xSPM struct
        [hReg,xSPM] = results_nogui('pv','spm',statdir,'con',spc{1},'thr',spc{2},'ptype',spc{3},'ext',spc{4},'disp',0);
        
        % Create Marsbar ROIs
        roipath = tgtdir;
        tmpstr1 = sprintf('_p%.3f',spc{2});
        tmpstr1(strfind(tmpstr1,'.')) = '_';
        if ischar(spc{3})
            tmpstr2 = spc{3};
        else
            tmpstr2 = num2str(spc{3});
        end
        
        rootn = ['SPM_',sprintf('con%.4d',spc{1}),tmpstr1,'_corr',tmpstr2,'_k',num2str(spc{4})];
        
        rootn = [xSPM.title,'_Ic',sprintf('%02.0f',xSPM.Ic),'_u',num2str(xSPM.u,3),'_k',num2str(xSPM.k),'_VOI'];
        rootn = strrep(rootn,' ','_');
        rootn = strrep(rootn,'.','_');
        rootn = strrep(rootn,':','_');
        
        
        mars_blobs2rois(xSPM, roipath, rootn);
        
    case 'image'
        
        switch outputtype
            case {'separate','cluster'}
                flags = 'i';        
            case 'mask'
                flags = 'c';
        end
        
        for ifls = 1:size(srcfls,1)
            tmpfl = deblank(srcfls(ifls,:));
            [p,rootn] = spm_fileparts(tmpfl);
            mars_img2rois(tmpfl,tgtdir,rootn,flags);
        end        
end
  
switch lower(varargin{1})
    case {'sphere','box'}

        voilist = [];
        for ivoi = 1:size(coord,1)
            voitag = [num2str(dims(end)),'-',num2str(coord(ivoi,1)),'_',num2str(coord(ivoi,2)),'_',num2str(coord(ivoi,3))];
            voilist = strvcat(voilist,spm_select('FPList',tgtdir,[voitag,'.*roi.mat$']));
        end
        
    case {'spm','image'}
        
        voilist = spm_select('FPList',tgtdir,['^',rootn,'.*roi.mat$']);
        
end

switch lower(varargin{1})
    case {'sphere','box','spm','image'}
        roispace = mars_space(spacedefimg);                                                        % Get space of preprocessed data
        
    switch lower(outputtype)
        case 'separate'

                for ivoi = 1:size(voilist,1)
                    [p n] = spm_fileparts(deblank(voilist(ivoi,:)));
                    voi = deblank(voilist(ivoi,:));
                    voi = maroi('load', voi);

                    imgfname = fullfile(tgtdir,[n,'.nii']);
                    save_as_image(voi, imgfname, roispace);
                end

        case {'cluster','mask'}

            switch lower(varargin{1})
                case {'sphere','box'}
                    clustfname = fullfile(tgtdir,[outputtype,'_',datestr(now,'yyyymmmdd'),'.nii']);                % Cluster map name
                case 'image'
                    clustfname = fullfile(tgtdir,[outputtype,'_',rootn,'_',datestr(now,'yyyymmmdd'),'.nii']);                % Cluster map name
                    clear roipath rootn tmpstr1 tmpstr2
                case 'spm'
                    fname = [xSPM.title,'_Ic',sprintf('%02.0f',xSPM.Ic),'_u',num2str(xSPM.u,3),'_k',num2str(xSPM.k)];
                    fname = strrep(fname,' ','_');
                    fname = strrep(fname,'.','_');
                    clustfname = fullfile(tgtdir,[fname,'_',datestr(now,'yyyymmmdd'),'.nii']); 
            end

            switch outputtype
                case 'cluster'
                    flags = 'i';        
                case 'mask'
                    flags = 'c';
            end

            mars_rois2img(voilist, clustfname, roispace, flags);

    end

end

switch lower(varargin{1})
    
    case 'lpba_mask'
        
        lpba_label = '/data/bram/ATLASES/LPBA40/LPBA40.SPM5.nifti/maxprob/lpba40.spm5.avg152T1.label.nii';
        
        for ifls = 1:size(srcfls,1)
            tmpfl = deblank(srcfls(ifls,:));
            [p,rootn,ext] = spm_fileparts(tmpfl);
            rtmpfl = bztbx_reslicemat(tmpfl,tmpfl,0);
            rlpba_label = bztbx_reslicemat(tmpfl,lpba_label,0);
            
            labelmask = rtmpfl.*rlpba_label;
            V = spm_vol(tmpfl);
            
            Vo = struct('fname',	fullfile(tgtdir,[rootn,'LPBA_labeled',ext]),...
                        'dim',		V.dim,...
                        'dt',		V.dt,...
                        'mat',		V.mat,...
                        'descrip', [V.descrip, ' - LPBA labeled']);
                       
            Vo = spm_write_vol(Vo,labelmask);
            
            bztbx_createmask('image',fullfile(tgtdir,[rootn,'LPBA_labeled',ext]),tmpfl,tgtdir,'cluster');
            delete(fullfile(tgtdir,[rootn,'LPBA_labeled',ext]));
            
            clear V Vo
            
            if k > 0
               
                % Load clustermap
                mapfl = spm_select('FPList',tgtdir,['^cluster_',rootn,'LPBA_labeled.*nii$']);
                V = spm_vol(mapfl);
                [Y,XYZ] = spm_read_vols(V);
                
                % Load clusterlabel file
                labelfl = spm_select('FPList',tgtdir,['^cluster_',rootn,'LPBA_labeled.*labels.mat$']);
                load(labelfl);
                
                
                survclust = [];
                
                for iclust = 1:max(Y(:))
                    if numel(find(Y(:) == iclust)) < k      % If current cluster smaller than cluster extent threshold
                       Y(find(Y(:) == iclust)) = 0;         % Remove cluster
                       eval(['!rm ',fullfile(tgtdir,[roi_info(iclust).label,'_roi.mat'])]);  % Remove cluster mat file
                    else
                        survclust = horzcat(survclust,iclust);
                    end
                end
                
                roi_info_new = [];
                
                Y_new = zeros(size(Y));
                for isc = 1:length(survclust)
                    Y_new(find(Y(:) == survclust(isc))) = isc;
                    roi_info_new(isc).label = roi_info(survclust(isc)).label;
                    roi_info_new(isc).number = isc;
                end
                
                clear roi_info
                roi_info = roi_info_new;
                save(fullfile(tgtdir,['cluster_',rootn,'LPBA_labeled_k',num2str(k),'_',datestr(now,'yyyymmmdd'),'_labels.mat']),'roi_info');
                
%                 K = find(Y(:) > 0);                         % Find cluster voxels
%                 XYZ(4,:) = 1;
%                 xyz = V.mat\XYZ(:,K);                       % Find voxels indices corresponding with cluster voxels
%                 xyz = xyz([1:3],:);
%                 A = spm_clusters(xyz);
%                 
%                 Y_new = zeros(size(Y));
%                 Y_new(K) = A;
                
                Vo = struct('fname',	fullfile(tgtdir,['cluster_',rootn,'LPBA_labeled_k',num2str(k),'_',datestr(now,'yyyymmmdd'),ext]),...
                        'dim',		V.dim,...
                        'dt',		V.dt,...
                        'mat',		V.mat,...
                        'descrip', [V.descrip, ' - LPBA labeled and cluster extent of ',num2str(k)]);
                       
                Vo = spm_write_vol(Vo,Y_new);
                
                
                % Remove original clustermap and clusterlabel file
                eval(['!rm ',mapfl]);
                eval(['!rm ',labelfl]);
            end
            clear tmpfl p rootn ext rtmpfl rlpba_label labelmask V Vo
        end
        
        
        
end






%
%         
%         
%         
%         
% for ivoi = 1:size(coord,1)
%         
%     switch char(fieldnames(voi(ivoi).voidef))
%         
%         case 'voisphere'
%             c = voi(ivoi).mnicoords;                                            % VOI center
%             r = voi(ivoi).voidef.voisphere.voispherespec;                       % VOI radius
%             d = sprintf('%0.1fmm radius sphere at [%0.1f %0.1f %0.1f]',r,c);    % VOI description
%             l = sprintf('voi_%02.0f_sphere_%0.0f-%0.0f_%0.0f_%0.0f',ivoi,r,c);                  % VOI label
%             o = maroi_sphere(struct('centre',c,'radius',r));                    % Create VOI object
%         case 'voibox'
%             c = voi(ivoi).mnicoords;                                            % VOI center
%             w = voi(ivoi).voidef.voibox.voiboxspec;                             % VOI dimensions
%             d = sprintf('[%0.1f %0.1f %0.1f] box at [%0.1f %0.1f %0.1f]',w,c);  % VOI description
%             l = sprintf('voi_%02.0f_box_w-%0.0f_%0.0f_%0.0f-%0.0f_%0.0f_%0.0f',ivoi,w,c);       % VOI label
%             o = maroi_box(struct('centre',c,'widths',w));                       % Create VOI object
%             
%     end
%     
%     o = descrip(o,d);                                                           % Add VOI description to VOI object
%     o = label(o,l);                                                             % Add VOI label to VOI object
%  
%     fn = fullfile(tgtdir,[label(o),'_roi.mat']);
%     saveroi(o,fn);                                                              % Save VOI
%     
% end
% 
% % Save each cluster/VOI as in a cluster image
% roi_list = [];
% for ivoi = 1:size(voi,2)
%     voitag = [num2str(voi(ivoi).mnicoords(1)),'_',num2str(voi(ivoi).mnicoords(2)),'_',num2str(voi(ivoi).mnicoords(3))];
%     roi_list = strvcat(roi_list,spm_select('FPList',tgtdir,[voitag,'.*roi.mat$']));
% end
% 
% img_name = fullfile(tgtdir,['clustermap_',datestr(now,'yyyymmmdd'),'.nii']);                % Cluster map name
% roi_space = mars_space(spacedefimg);                                                        % Get space of preprocessed data        
% flags = 'i';                                                                                % Clustermap
% mars_rois2img(roi_list, img_name, roi_space, flags);
% 
% 
% 
% % Save each cluster/VOI as separate image
% roi_space = mars_space(spacedefimg);                                                        % Get space of preprocessed data        
% 
% for iroi = 1:size(roi_list,1)
%     [p n e v] = spm_fileparts(deblank(roi_list(iroi,:)));
%     roi = deblank(roi_list(iroi,:));
%     img_name = fullfile(tgtdir,[n,'.nii']);
%     save_as_image(roi, img_name, roi_space);
% end
% 
% 
% 
% 
