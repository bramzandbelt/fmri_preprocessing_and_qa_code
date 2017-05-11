function bztbx_bpfilt(varargin)
% ________________________________________________________________________
% BAND-PASS FILTER Performs band-pass filtering for resting-state fMRI data
%
% FORMAT: bztbx_bpfilt(files,isi,lowpb,highpb,norder);
% DESCRIPTION       Standard band-pass filtering
% files             - filenames of realigned scans (char array)
% isi               - interscan interval (scalar, in seconds)
% lowpb             - lowest pass band (scalar, in Hertz)
% highpb            - highest pass band (scalar, in Hertz)
% norder            - filter order
%
% FORMAT: bztbx_bpfilt(files,isi,lowpb,highpb,norder,rpfile);
% DESCRIPTION       Same as standard, but also truncates realignment
%                   parameter file for 1st level fMRI analysis
% rpfile            - realignment parameter file (char array)
% ________________________________________________________________________
% Bram B. Zandbelt, August 2008
% b.b.zandbelt@umcutrecht.nl
%
% version 1.0

if nargin == 5
    ts = varargin{1};
    isi = varargin{2};
    lowpb = varargin{3};
    highpb = varargin{4};
    norder = varargin{5};
elseif nargin == 6
    ts = varargin{1};
    isi = varargin{2};
    lowpb = varargin{3};
    highpb = varargin{4};
    norder = varargin{5};
    rpfile = varargin{6};
else
    error('Wrong number of input arguments. Type ''help bztbx_bpfilt'' for further assistance.');
end
   
if isempty(ts)
    [ts] = spm_select(inf,'image','Select Nifti files');
end

% Loading data
V = spm_vol(ts);
nScan = size(V,1);
xdim = V(1).dim(1);
ydim = V(1).dim(2);
zdim = V(1).dim(3);

Mask = ones(xdim,ydim,zdim);

% Set-up high pass filter
HPF.HParam = nScan*isi/2+1;                                 % Period of half the experiment duration + 1
HPF.row = 1:nScan;
HPF.RT = isi;
HPF = spm_filter(HPF);

% Set-up band-pass filter
Fs=1/isi;
[filt1,filt2]=fir1(norder,[lowpb/Fs/2,highpb/Fs/2],'bandpass');


% Get plane coordinates
xords = (1:xdim)'*ones(1,ydim); xords = xords(:)';
yords = ones(xdim,1)*(1:ydim);  yords = yords(:)';

% Initialise output files
Vbpf(1:nScan) = deal(struct(...
    'fname',    [],...
    'dim',      V(1).dim,...
    'dt',       V(1).dt,...
    'mat',      V(1).mat,...
    'pinfo',    V(1).pinfo,...
    'descrip',  ''));
for i = 1:nScan
    [p n e v] = spm_fileparts(V(i).fname);
    Vbpf(i).fname = fullfile(p,['f',n,e,v]);
    Vbpf(i).descrip   = [V(i).descrip,' - band-pass (',num2str(lowpb),' - ',num2str(highpb),' Hz) filtered'];
    clear p n e v      
end

Vbpf = spm_create_vol(Vbpf);

fprintf('%-40s: %30s','Band-pass filtering','...computing')    %-#
spm_progress_bar('Init',100,'Band-pass filtering','');

% CODE BELOW WILL BE USED ELSEWHERE
% % For plotting
% tsrawsl = nan(nScan,zdim);  % raw timeseries
% tshpfsl = nan(nScan,zdim);  % high-pass filtered time series
% tsbpfsl = nan(nScan,zdim);  % band-pass filtered time series
% tsfilsl = nan(nScan,zdim);  % filtered and mean added time series
% 
for z = 1:zdim
    
    % current plane-specific parameters
    zords   = z*ones(xdim*ydim,1)';             % plane Z coordinates
    Cm = reshape(Mask(:,:,z),[1,xdim*ydim]);    
    Cm = find(Cm>0);  
    
    % print progress information in command window
    str   = sprintf('Plane %3d/%-3d',z,zdim);
    fprintf('\r%-40s: %30s',str,' ')

    % construct list of voxels in this block
    xyz   = [xords; yords; zords];              % voxel coordinates
    nVox  = size(xyz,2);                        % number of voxels in plane

    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...reading data')  %-#
    
    % Get the data
    Y     = zeros(nScan,nVox);
    for i = 1:nScan
        Y(i,Cm)  = spm_get_data(V(i).fname,xyz(:,Cm));
    end
    
    tsrawsl(:,z) = nanmean(Y(:,Cm),2);
    
    b0 = nanmean(Y,1);
    
    Y = spm_filter(HPF,(Y - kron(ones(nScan,1),b0)));       % High-pass filter mean removed ts
%     tshpfsl(:,z) = nanmean(Y(:,Cm),2);
    
    
    Y = filter(filt1,filt2,Y);                              % Band-pass filter data
%     tsbpfsl(:,z) = nanmean(Y(:,Cm),2);
   
    Y = Y + kron(ones(nScan,1),b0);                         % Adding mean again
%     tsfilsl(:,z) = nanmean(Y(:,Cm),2);
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...saving plane')  %-#
 
    % Write filtered data
    for i = 1:nScan
        jj   = NaN*ones(xdim,ydim);
        jj(:) = Y(i,:);
        spm_write_plane(Vbpf(i), jj, z);
    end

    % Report progress
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done')
    spm_progress_bar('Set',100*(z/zdim));
    
end
 
% tsraw = nansum(tsrawsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average raw time series across all mask voxels
% tshpf = nansum(tshpfsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average raw time series across all mask voxels
% tsbpf = nansum(tsbpfsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average raw time series across all mask voxels
% tsfil = nansum(tsfilsl.*kron(ones(nScan,1),(nMVoxSl/nMVox)),2); % Average raw time series across all mask voxels
% 
% if nargin == 6
%     [p n e v] = spm_fileparts(rpfile);
%     rp = load(rpfile);
%     rp_trunc = rp(1+norder*2:end,:);
%     save(fullfile(p,['t',n,e,v]),'rp_trunc','-ASCII');
% end 
%     
% % Generating figure
% fig = spm_figure('FindWin','Graphics');
% if isempty(fig)
%     fig = spm_figure('Create','Graphics');
% end;
% set(0,'CurrentFigure',fig);
% spm_figure('Clear','Graphics');
% 
% t = ((1:nScan-2*norder)-1)/Fs;                              % Time in seconds
% 
% subplot(4,2,1)
% plot(t,tsraw(2*norder+1:nScan));
% title('Signal - Raw');
% 
% subplot(4,2,2)
% [Pxx_raw,f_raw] = pwelch(tsraw(2*norder+1:nScan),[],[],[],Fs);
% plot(f_raw,Pxx_raw);
% title('Spectrum - Raw');
% 
% subplot(4,2,3)
% plot(t,tshpf(2*norder+1:nScan));
% title('Signal - HPF');
% 
% subplot(4,2,4)
% [Pxx_hpf,f_hpf] = pwelch(tshpf(2*norder+1:nScan),[],[],[],Fs);
% plot(f_hpf,Pxx_hpf);
% title('Spectrum - HPF');
% 
% subplot(4,2,5)
% plot(t,tsbpf(2*norder+1:nScan));
% title('Band pass filtered');
% 
% subplot(4,2,6)
% [Pxx_bpf,f_bpf] = pwelch(tsbpf(2*norder+1:nScan),[],[],[],Fs);
% plot(f_bpf,Pxx_bpf);
% title('Spectrum - BPF');
% 
% subplot(4,2,7)
% plot(t,tsfil(2*norder+1:nScan));
% title('Mean added');
% 
% subplot(4,2,8)
% [Pxx_ma,f_ma] = pwelch(tsfil(2*norder+1:nScan),[],[],[],Fs);
% plot(f_ma,Pxx_ma);
% title('Spectrum - Mean added');
