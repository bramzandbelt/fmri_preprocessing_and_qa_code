function bztbx_tsnr(varargin)
% _________________________________________________________________________
% TSNR Computes voxel-wise temporal signal-to-noise ratio
%
% FORMAT        bztbx_tsnr('basic',ts,isi,tgtdir,mask);
% DESCRIPTION   tsnr computation of raw data
% ts            - char array of filenames (or leave empty to select)
% isi           - interscan interval (in seconds)
% tgtdir        - char array of target directory path
% mask          - char array of filename of mask to which analysis should be 
%               restricted (optional)
%
% FORMAT        bztbx_tsnr('hpf',ts,isi,hp_cutoff,tgtdir,mask);
% DESCRIPTION   tsnr computation of high-pass filtered data
% hp_cutoff     - highpass cut-off period (in seconds)
% tgtdir        - char array of target directory path
%
% FORMAT:       bztbx_tsnr('writehpf',ts,isi,hp_cutoff,tgtdir,mask);
% DESCRIPTION   as 'hpf', but also writes filtered data
%
% _________________________________________________________________________
% Bram B. Zandbelt, July 2008
% b.b.zandbelt@umcutrecht.nl
%
% version 1.0

warning off % To disable 'Divided by zero' message

% Check inputs
if nargin == 0
    error('No input arguments provided. Type help tsnr for usage')
end

switch lower(varargin{1}) 
    case 'basic'
        ts = varargin{2};
        isi = varargin{3};
        tgtdir = varargin{4};
        if nargin == 5;
            mask = varargin{5};
        elseif nargin > 5
            error('Wrong number of input arguments. Type ''help bztbx_tsnr'' for usage')
        end
    case {'hpf','writehpf'}
        ts = varargin{2};
        isi = varargin{3};
        hp_cutoff = varargin{4};
        tgtdir = varargin{5};
        if nargin == 6;
            mask = varargin{6};
        elseif nargin > 6
            error('Wrong number of input arguments. Type ''help bztbx_tsnr'' for usage')
        end
end

if isempty(ts)
    [ts] = spm_select(inf,'image','Select Nifti files');
end

if isempty(tgtdir)
    [tgtdir] = spm_select(1,'dir','Select target directory');
end

V = spm_vol(ts);
nScan = size(V,1);
xdim = V(1).dim(1);
ydim = V(1).dim(2);
zdim = V(1).dim(3);

% Create or load brain mask based on first fMRI volume
if exist('mask','var') == 1
    VM = spm_vol(mask);
    Mask = spm_read_vols(VM);
    nMVox = numel(find(Mask(:) > 0));
    nMVoxSl = zeros(1,zdim);
    for i = 1:zdim
        nMVoxSl(1,i) = numel(find(Mask(:,:,i) > 0));
    end
else
    Mask = bztbx_automask('fast',V(1).fname);
    nMVox = numel(find(Mask(:) > 0));
    nMVoxSl = zeros(1,zdim);
    for i = 1:zdim
        nMVoxSl(1,i) = numel(find(Mask(:,:,i) > 0));
    end
end

% Mask = bztbx_automask(V(1).fname);
% [M,Mfile] = bztbx_automask('robust',srcimg,spacedefimg);
% nMVox = numel(find(Mask(:) > 0));
% nMVoxSl = zeros(1,zdim);
% for i = 1:zdim
%     nMVoxSl(1,i) = numel(find(Mask(:,:,i) > 0));
% end


switch lower(varargin{1})
    case {'hpf','writehpf'}
        % Set-up high pass filter
        HPF.HParam = hp_cutoff;
        HPF.row = 1:nScan;
        HPF.RT = isi;
        HPF = spm_filter(HPF);
end

% Get plane coordinates
xords = (1:xdim)'*ones(1,ydim); xords = xords(:)';
yords = ones(xdim,1)*(1:ydim);  yords = yords(:)';

[p n e v] = spm_fileparts(V(1).fname);

if isempty(tgtdir)
    tgtdir = fullfile(p,'QAfiles/');
end

if exist(tgtdir,'dir') ~= 7
    mkdir (tgtdir)
end

% Initialise output files
Vtmean = struct(...
        'fname',    fullfile(tgtdir,['tmean_',n,e,v]),...
        'dim',      V(1).dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      V(1).mat,...
        'pinfo',    [1 0 0]',...
        'descrip',  [V(1).descrip,' - temporal mean signal map']);
Vtmean = spm_create_vol(Vtmean);
    
Vtsd = struct(...
        'fname',    fullfile(tgtdir,['tsd_',n,e,v]),...
        'dim',      V(1).dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      V(1).mat,...
        'pinfo',    [1 0 0]',...
        'descrip',  [V(1).descrip,' - temporal standard deviation signal map']);
Vtsd = spm_create_vol(Vtsd);    

Vtsnr = struct(...
        'fname',    fullfile(tgtdir,['tsnr_',n,e,v]),...
        'dim',      V(1).dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      V(1).mat,...
        'pinfo',    [1 0 0]',...
        'descrip',  [V(1).descrip,' - temporal signal-to-noise ratio map']);
Vtsnr = spm_create_vol(Vtsnr);  

clear p n e v

switch lower(varargin{1})
    case {'writehpf'}
     Vhpf(1:nScan) = deal(struct(...
        'fname',    [],...
        'dim',      V(1).dim,...
        'dt',       V(1).dt,...
        'mat',      V(1).mat,...
        'pinfo',    V(1).pinfo,...
        'descrip',  ''));
    for i = 1:nScan
        [p n e v] = spm_fileparts(V(i).fname);
        Vhpf(i).fname = fullfile(p,['hpf',num2str(hp_cutoff),n,e,v]);
        Vhpf(i).descrip   = [V(i).descrip,' - high-pass (',num2str(hp_cutoff),'s) filtered'];
        clear p n e v      
    end

    Vhpf = spm_create_vol(Vhpf);
    
end


% Determine grand mean scaling factor
for i = 1:nScan
    g(i) = spm_global(V(i));
end
gSF   = 100./g';


fprintf('%-40s: %30s','Time series descriptives','...computing')    %-#

str = ['tsnr: ',varargin{1}];
spm_progress_bar('Init',100,str,'');

tsrawsl = nan(nScan,zdim);
tshpfsl = nan(nScan,zdim);

for z = 1:zdim
    
    % current plane-specific parameters
    zords   = z*ones(xdim*ydim,1)'; %-plane Z coordinates
    Cm = reshape(Mask(:,:,z),[1,xdim*ydim]);    
    Cm = find(Cm>0);                            % Mask indices
    
    % print progress information in command window
    str   = sprintf('Plane %3d/%-3d',z,zdim);
    fprintf('\r%-40s: %30s',str,' ')

    % construct list of voxels in this block
    xyz   = [xords; yords; zords];              %-voxel coordinates
    nVox  = size(xyz,2);                        %-number of voxels in plane

    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...reading data')  %-#
    
    % Get the data
    Y     = zeros(nScan,nVox);
    for i = 1:nScan
        Y(i,Cm)  = spm_get_data(V(i).fname,xyz(:,Cm));
    end

    
    tsrawsl(:,z) = nanmean(Y(:,Cm),2);
       
    % High-pass filter data
    switch lower(varargin{1})
        case {'hpf','writehpf'}
            Y    = spm_filter(HPF,Y.*kron(gSF,ones(1,nVox)));          % Data is GM scaled before HPF
            tshpfsl(:,z) = nanmean(Y(:,Cm),2);
        otherwise
            HPF.HParam = 128;
            HPF.row = 1:nScan;
            HPF.RT = isi;
            HPF = spm_filter(HPF);
            Y    = spm_filter(HPF,Y.*kron(gSF,ones(1,nVox)));          % Data is GM scaled before HPF
            tshpfsl(:,z) = nanmean(Y(:,Cm),2);
%             Yhpf = spm_filter(HPF,Y.*gSF);                        % OLD CODE
%             tshpfsl(:,z) = nanmean(Yhpf(:,Cm),2);
%             clear Yhpf
    end
   
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...saving plane')  %-#
 
    switch lower(varargin{1})
        case 'writehpf'
   
        % Write filtered data
        for i = 1:nScan
            jj   = NaN*ones(xdim,ydim);
            jj(:) = Y(i,:);
            spm_write_plane(Vhpf(i), jj, z);
        end
        
    end

    % Write temporal mean map
    jj   = NaN*ones(xdim,ydim);
    jj(:) = nanmean(Y,1);
    spm_write_plane(Vtmean,jj,z);

    % Write temporal standard deviation map
    jj   = NaN*ones(xdim,ydim);
    jj(:) = nanstd(Y,1);
    spm_write_plane(Vtsd,jj,z);

    % Write temporal signal-to-noise ratio map
    jj   = NaN*ones(xdim,ydim);
    jj(:) = nanmean(Y,1)./nanstd(Y,1);
    spm_write_plane(Vtsnr,jj,z);
    
    % Report progress
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done')
    spm_progress_bar('Set',100*(z/zdim));

end

% Write mean raw and hpf signal per slice per volume
tsraw = nansum(tsrawsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average raw time series across all mask voxels
tshpf = nansum(tshpfsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average hpf time series across all mask voxels

tspsc = abs(tshpf*100/nanmean(tshpf)-100);                          % Global percent signal change
tspscsl = tshpfsl.*kron(ones(nScan,1),100./nanmean(tshpfsl,1)) - 100; % Local (i.e. slice) percent signal change

[p n] = spm_fileparts(V(1).fname);
save(fullfile(tgtdir,['tsrawsl_',n,'.mat']),'tsrawsl');
save(fullfile(tgtdir,['tshpfsl_',n,'.mat']),'tshpfsl');
save(fullfile(tgtdir,['tsraw_',n,'.mat']),'tsraw');
save(fullfile(tgtdir,['tshpf_',n,'.mat']),'tshpf');
save(fullfile(tgtdir,['tspsc_',n,'.mat']),'tspsc');
save(fullfile(tgtdir,['tspscsl_',n,'.mat']),'tspscsl');