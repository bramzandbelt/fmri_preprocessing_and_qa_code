function  varargout = bztbx_automask(varargin)
% _________________________________________________________________________
% AUTOMASK Computes a brain mask on the basis of a single fMRI volume
%
% There are two methods, one based on Paul Mazaika's ArtRepair Toolbox 
% (Mazaika et al. (2005), Hum Brain Mapp), the other on Stephen Smith's 
% Brain Extraction Tool (Smith (2002) Hum Brain Mapp). The ArtRepair method
% is fast (and therefore referred to as 'fast') and performs pretty well on 
% fMRI data, but may includes the eyes and parts of the skull. The Brain 
% Extraction Tool method, provided that a T1-image is used, is more robust 
% (and therefore referred to as 'robust'), but also time-consuming. 
%
%
% FORMAT        [M,Mfile] = bztbx_automask('fast',srcimg);
% DESCRIPTION   Computes brain mask with defaults settings using
%               ArtRepair's method
% srcimg        - char array of filename
% M             - mask
% Mfile         - mask file name
%
% FORMAT        [M,Mfile] = bztbx_automask('fast',srcimg,threshold,writevol);
% DESCRIPTION   Computes brain mask with user-specifed settings using
%               ArtRepair's method
% threshold     - threshold value for calculating the mask
% writevol      - write (1) or do not write auto mask file
%
% FORMAT        [M,Mfile] = bztbx_automask('robust',srcimg,spacedefimg);
% DESCRIPTION   Computes brain mask using the Brain Extraction Tool method
% - srcimg      - char array of filename
% - spacedefimg - char array of filename of space defining image (optional)
%
% _________________________________________________________________________
% Bram Zandbelt, October 2009
% b.b.zandbelt@umcutrecht.nl
%
% version 1.1
%


switch varargin{1}
    case 'fast'

        % Code is based on ArtRepair Toolbox
        FMR = 0.015;   % False Mask Percentage. 

        if nargin == 1
            srcimg = spm_select(1,'image','Select image as source for automask.',[],pwd);
            Threshold = -1;  %  Adaptive mask by default
            WriteVol = 1;    %  Writes the mask out for inspection.
        elseif nargin == 2
            Threshold = -1;  %  Adaptive mask by default.
            WriteVol = 1;    %  Writes the mask out for inspection.
            srcimg = varargin{2};
        elseif nargin == 4
            srcimg = varargin{2};
            Threshold = varargin{3};
            WriteVol = varargin{4};
        else
            error('Wrong number of input arguments provided. Type ''help bztbx_automask'' for usage')
        end


        % Get the image data.
            V = spm_vol(srcimg(:,:));  % Input could be a matrix; only need one image.
            Y = spm_read_vols(V);
        % Fill in the small holes and reduce the noise spikes.    
            Y = smooth3(Y);  % default 3x3x3 box smoothing.
            Yr = max(max(max(Y))) - min(min(min(Y)));  % previously used range.

            % User defined mask threshold
        if Threshold > 0  % Make the mask directly
            % Array temp is logical array with 1's and 0's
            temp(:,:,:) = (Y(:,:,:)>Yr*Threshold);
        end

        % Adaptive Mask Threshold
        if ( Threshold == -1 )   % Find a threshold that discards three outer faces.
            % Use gray matter density as lower limit- count is 400.
            Tlow = fix(0.2*Yr);  Thigh = fix(0.4*Yr); Tskip = max(fix(Tlow/20),1);
            for Tbar = Tlow:Tskip:Thigh   % 400:20:800
                temp(:,:,:) = (Y(:,:,:) > Tbar);
            % Count the number of mask points in the far faces of the volume
                 xdim = size(Y);
                 count1 = sum(sum(temp(:,:,1)));
                 count2 = sum(sum(temp(:,:,xdim(3))));
                 count3 = sum(sum(temp(:,1,:)));
                 count4 = sum(sum(temp(:,xdim(2),:)));
                 count5 = sum(sum(temp(1,:,:)));
                 count6 = sum(sum(temp(xdim(1),:,:)));
             % Always have one face with large counts, sometimes have 2 such faces.
                 Xbar = [ count1 count2 count3 count4 count5 count6 ];
                 Ybar = sort(Xbar);
                 countC = Ybar(1) + Ybar(2);  % the two smallest face counts
                 % Number of voxels on 3 faces is approximately:
                 nvox = xdim(1)*xdim(2) + xdim(2)*xdim(3) + xdim(1)*xdim(3);
                 if ( countC < FMR*nvox )  
                     break;   % Exit the For loop, current Tbar is good enough.
                 end

             end

        end
        if Threshold == -1
            if Tbar >= Thigh-Tskip
                disp('Automask program failed.   Try choosing a mean image,')
                disp(' or manually set a threshold. Type ''help bztbx_automask'' for further assistance.')
                return
            end
        end


        if ( WriteVol == 1 )     
            v = V;  % preserves the header structure
            [dirname, xname, xext ] = fileparts(V.fname);
            artifname = ['automask_' xname xext];
            automaskfile = fullfile(dirname,artifname);
            v.fname = automaskfile;
            noscale_write_vol(v,Y);  
        end

        % Outputs
        varargout{1} = temp;
        varargout{2} = automaskfile;
        
    case 'robust'
        
        if nargin == 1
            srcimg = spm_select(1,'image','Select image as source for automask.',[],pwd);
        elseif nargin == 2
            srcimg = varargin{2};
            opts = '-R -m';
        elseif nargin == 3
            srcimg = varargin{2};
            spacedefimg = varargin{3};
            opts = '-R -m';
        else
            error('Wrong number of input arguments provided. Type ''help bztbx_automask'' for usage')
        end
        
        try
            delete(fullfile(p,[n,'_brain_mask.nii']));
            delete(fullfile(p,[n,'_brain.nii']));
        catch
        end
        
        [p n] = spm_fileparts(srcimg);
        cmd = ['/usr/local/fsl/bin/bet ',fullfile(p,n),' ',fullfile(p,[n,'_brain']),' ',opts];
        system(cmd); % Run Brain Extraction Tool
    
        if exist(fullfile(p,'masknz.nii.gz'),'file') == 2
            delete(fullfile(p,'masknz.nii.gz'));
        end
        
        gunzip(cellstr(spm_select('FPList',p,[n,'.*brain.*gz$']))); % Unzip the BET files
        delete(fullfile(p,[n,'*brain*gz'])); % Remove zipped data
        
        if exist('spacedefimg','var') == 1
            P = strvcat(spacedefimg,spm_select('FPList',p,['^',n,'.*brain_mask.nii$'])); % Identfy files
            spm_reslice(P,{0 0 0 1 [0 0 0]}); % Reslice the brain mask to dimensions of space defining image
            
            % Outputs
            varargout{1} = spm_read_vols(spm_vol(spm_select('FPList',p,['^r',n,'.*brain_mask.nii$'])));
            varargout{2} = spm_select('FPList',p,['^r',n,'.*brain_mask.nii$']);
        else
            % Outputs
            varargout{1} = spm_read_vols(spm_vol(spm_select('FPList',p,['^',n,'.*brain_mask.nii$'])));
            varargout{2} = spm_select('FPList',p,['^',n,'.*brain_mask.nii$']);
        end
        clear p n e v gzfile cmd
        
end
         
%---------------------------------------------------------------
% Create and write image without the scale and offset steps  
% This function is spm_write_vol without error checking and scaling.
function noscale_write_vol(V,Y)
V = spm_create_vol(V);
for p=1:V.dim(3),
    V = spm_write_plane(V,Y(:,:,p),p);
end;
%V = spm_close_vol(V);  % not for SPM5  