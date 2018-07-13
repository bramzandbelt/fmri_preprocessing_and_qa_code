function bztbx_extractts(varargin)
% ________________________________________________________________________
% TIME SERIES EXTRACTION 
%
% FORMAT: bztbx_extractts('raw',ts,mask,tgtdir,tag);
% DESCRIPTION       Extracts raw time series
% ts                - filenames of fMRI time series (char array)
% mask              - filenames of masks to extract data from (char array)
% tgtdir            - filename of target directory (char array)
% tag               - suffix (char array)
%
% FORMAT: bztbx_extractts('detrended',ts,mask,tgtdir,tag);
% DESCRIPTION       Extracts detrended (mean = 0) time series.
%
% FORMAT: bztbx_extractts('standardized',ts,mask,tgtdir,tag);
% DESCRIPTION       Extracts standardized (mean = 0, std = 1) time series.
%
% FORMAT: bztbx_extractts('hpf',ts,mask,tgtdir,tag,isi,hpf);
% DESCRIPTION       Extracts grand-mean scaled, high-pass filtered time
%                   series
% isi               - interscan interval, in seconds (scalar)
% hpf               - high-pass filter cut-off, in seconds (scalar)
%
% ________________________________________________________________________
% Bram B. Zandbelt, October 2009
% b.b.zandbelt@umcutrecht.nl
%
% version 1.0

% Check inputs
if nargin == 0
    error('No input arguments provided. Type ''help bztbx_extractts'' for usage')
end

switch lower(varargin{1}) 
    case {'raw','detrended','standardized'}
        ts = varargin{2};
        mask = varargin{3};
        tgtdir = varargin{4};
        tag = varargin{5};
    case {'hpf'}
        ts = varargin{2};
        mask = varargin{3};
        tgtdir = varargin{4};
        tag = varargin{5};
        isi = varargin{6};
        hp_cutoff = varargin{7};
end

if isempty(ts)
    ts = spm_select(inf,'image','Select Nifti files');
end

if isempty(mask)
    mask = spm_select(inf,'image','Select mask image(s)');
end

if isempty(tgtdir)
    tgtdir = spm_select(1,'dir','Select target directory');
end

% Load ts data
V = spm_vol(ts);
nScan = size(V,1);
xdim = V(1).dim(1);
ydim = V(1).dim(2);
zdim = V(1).dim(3);

% Load mask data
VM = spm_vol(mask);
tsmean = zeros(nScan,size(VM,1));

% Progress bar setup
spm_progress_bar('Init',100,'Time series extraction','');


for iMask = 1:size(VM,1)

    Mask = spm_read_vols(VM(iMask));
    nMVox = numel(find(Mask(:) > 0));
    nMVoxSl = zeros(1,zdim);
    for i = 1:zdim
        nMVoxSl(1,i) = numel(find(Mask(:,:,i) > 0));        % For weighted average
    end

    switch lower(varargin{1})
        case 'hpf'

        % Set-up high pass filter
        HPF.HParam = hp_cutoff;
        HPF.row = 1:nScan;
        HPF.RT = isi;
        HPF = spm_filter(HPF);

        % Determine grand mean scaling factor
        for i = 1:nScan
            g(i) = spm_global(V(i));
        end
        gSF   = 100./g';

    end

    % Get plane coordinates
    xords = (1:xdim)'*ones(1,ydim); xords = xords(:)';
    yords = ones(xdim,1)*(1:ydim);  yords = yords(:)';

    [p n e v] = spm_fileparts(V(1).fname);
    
    countz = 0;
    for z = find(nMVoxSl > 0)
        
        countNz = numel(find(nMVoxSl > 0));
        countz = countz + 1;
        
        % current plane-specific parameters
        zords   = z*ones(xdim*ydim,1)';             % Plane Z coordinates
        Cm = reshape(Mask(:,:,z),[1,xdim*ydim]);    
        Cm = find(Cm>0);                            % Mask indices

        % construct list of voxels in this block
        xyz   = [xords; yords; zords];              % Voxel coordinates
        nVox  = size(xyz,2);                        % Number of voxels in plane

        % Get the data
        Y     = zeros(nScan,nVox);
        for i = 1:nScan
            Y(i,Cm)  = spm_get_data(V(i).fname,xyz(:,Cm));
        end

        % GM scale and HP filter if required
        switch lower(varargin{1})
            case 'hpf'
                Y = spm_filter(HPF,Y.*kron(gSF,ones(1,nVox)));
        end
        
        tsmean(:,iMask) = nansum([tsmean(:,iMask), nMVoxSl(z)/nMVox*nanmean(Y(:,Cm),2)],2);

        spm_progress_bar('Set',100*((iMask - 1)/size(VM,1) + 1/size(VM,1)*countz/countNz));
        
    end

    % Detrend or standardize if required
    switch lower(varargin{1}) 
        case 'detrended'
            tsmean(:,iMask) = spm_detrend(tsmean(:,iMask));
        case 'standardized'
            tsmean(:,iMask) = spm_detrend(tsmean(:,iMask))./std(tsmean(:,iMask),0,1);
    end
    
%     [p n] = spm_fileparts(V(1).fname);
%     save(fullfile(tgtdir,['mean_',lower(varargin{1}),'_ts_',n,'.mat']),'tsmean'); 
    
end

[p n] = spm_fileparts(V(1).fname);
save(fullfile(tgtdir,['mean_',lower(varargin{1}),'_ts_',n(1:end-5),'_',tag,'.txt']),'tsmean','-ascii');