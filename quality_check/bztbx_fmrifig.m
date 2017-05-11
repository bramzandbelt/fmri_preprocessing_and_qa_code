function varargout = bztbx_fmrifig(varargin)
% ________________________________________________________________________
% GENERATES ACTIVATION/CLUSTER MAPS OVERLAID ON ANATOMICAL SCANS
%
% FORMAT bztbx_fmrifig(t1img,overlay1,...,settings)
% DESCRIPTION  
% t1img         1x4 cell array, containing
%   type        - overlay type:'truecolor','split','contour' (default='truecolor')
%   fname       - filename (obligatory)
%   range       - range (default={mn mx})
%   cmap        - colormap name (default='gray')
% 
% overlay1      1x4 cell array (optional), contents depend on type:
%   type        - 'split' (SPMs)
%   fname       - filename of SPM.mat file
%   statspec    - cell array specifying statistical analysis
%       con       - contrast number
% 	    p         - significance level (e.g. .05)
%       ptype     - multiple comparison type ('none','FDR','FWE')
% 	    k         - cluster extent threshold
%       mspm      - path to the SPM file of the mask
%       mcon      - contrast number for masking (optional)
%       mp        - contrast threshold for masking (optional, MUST be p-value)
%       mptype    - inclusive (0) or exclusive (1)
%                  e.g. {1,.001,'none'15} or {1,.005,'none'39,1,.005,0}               
%
%   cmap        - colormap name (default:'autumn')
%
%   type        - 'split' (imgages, e.g. masks, cluster maps)
%   fname       - filename of image
%   range       - range (default:{mn mx})
%   cmap        - colormap name (default:'actc')
%
%   type        - 'truecolor' (transparent images)
%   fname       - filename of image
%   range       - range (default:{mn mx})
%   cmap        - colormap name (default:'actc')
%   prop        - proportion of intensity (optional, use negative values (e.g. -.25) in combination with a gray colormap with range [-1 0] to
%                 create a transparency mask)
%
%   type        - 'contour'
%   fname       - filename of image (mask)
%   linespec    - specification of the line (e.g. 'w-', see MATLAB's linespec)
%
% overlayN      optional, see overlay1
% 
% settings      1x9 cell array, containing
%   orient      - orientation: 'axial','coronal','sagittal' (default='axial')
%   slice       - array containing slices to display (in mm,default=[-40 -20 0 20 40])
%   cols        - number of columns in output (default=5)
%   dispcbar    - color bars to display (default=[])
%   lab         - labels to display (default = []): 'none',[]
%   res         - resolution (in dpi,default=300)
%   fname       - filename of output (default=fullfile(pwd,[fMRIfig_[date].tiff])
%   papwidth    - paperwidth (in cm,default=17)
%   visible     - yes (1) or no (0). if 0, the slover object is put into
%               workspace
%
% FORMAT settings = bztbx_fmrifig('defsettings')
% DESCRIPTION   Returns cell array with default settings
% ________________________________________________________________________
% Bram B. Zandbelt, January 2010
% b.b.zandbelt@umcutrecht.nl
%
% version 1.0

if ischar(varargin{1})
    switch lower(varargin{1})
        case 'defsettings'
            settings{1} = 'axial';
            settings{2} = [-20:8:68];
            settings{3} = 8;
            settings{4} = [];
            settings{5} = 'none';
            settings{6} = 450;
            settings{7} = [];
            settings{8} = 17;
            settings{9} = 1;
            varargout{1} = settings;
            return
    end
end


settings = varargin{nargin};

% addpath /mnt/home/bram/m-files/thomas_nichols/      % For CorrClusTh function

so = slover();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure T1-weighted image %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Layer type
if isempty(varargin{1}{1})
    so.img(1).type = 'truecolor';
else 
    so.img(1).type = varargin{1}{1};
end

% Header info of layer
if isempty(varargin{1}{2})
    mesg = ['Select image file for layer 1'];
    so.img(1).vol = spm_vol(spm_select(1,'image',mesg));
else
    so.img(1).vol = spm_vol(varargin{1}{2});
end

% Range 
if isempty(varargin{1}{3})
    [mx mn] = slover('volmaxmin', so.img(1).vol);
    so.img(1).range = [mn mx];
elseif isempty(varargin{1}{3}{1}) | isempty(varargin{1}{3}{2})
            [mx mn] = slover('volmaxmin', so.img(1).vol);
            rng = [];
            if isempty(varargin{1}{3}{1})
                rng(1) = mn;
            else
                rng(1) = varargin{1}{3}{1};
            end
            if isempty(varargin{1}{3}{2})
                rng(2) = mx;
            else
                rng(2) = varargin{1}{3}{2};
            end
            so.img(1).range = rng;   
else
    so.img(1).range = [varargin{1}{3}{1} varargin{1}{3}{2}];
end

% Color map
if isempty(varargin{1}{4})
    so.img(1).cmap = slover('getcmap','gray'); % default
else
    if ischar(varargin{1}{4})
        so.img(1).cmap = slover('getcmap',varargin{1}{4});
    elseif isnumeric(varargin{1}{4})
        so.img(1).cmap = varargin{1}{4};
    end
end

% Other parameters
% so.img(1).range = [mn mx];
so.img(1).outofrange = {0 64};
% so.img(1).nancol = [0 0 0]; % For black background
so.img(1).nancol = [1 1 1]; % For white background

%%%%%%%%%%%%%%%%%%%%%%
% Configure overlays %
%%%%%%%%%%%%%%%%%%%%%%

for ilayer = 2:nargin-1
    
    switch lower(varargin{ilayer}{1})
        
        case 'split'
            % Layer type
            so.img(ilayer).type = 'split';
            
            % Determine whether file is image or mat file
            if isempty(varargin{ilayer}{2})
                mesg = ['Select image or SPM.mat file for layer 1'];
                fl = spm_select(1,'any',mesg);
            else
                fl = varargin{ilayer}{2};
            end
            [p n e] = spm_fileparts(fl);
            
            switch e
                case {'.img','.hdr','.nii'}
                    % Header info of layer
                    v_b = spm_vol(fl);
                    [y_b xyz_b] = spm_read_vols(v_b);
                    so.img(ilayer).vol = slover('matrix2vol', y_b, v_b.mat);
                    
                    % Range
                    if isempty(varargin{ilayer}{3})
                        [mx mn] = slover('volmaxmin', so.img(ilayer).vol);
                        so.img(ilayer).range = [mn mx];
                    elseif isempty(varargin{ilayer}{3}{1}) | isempty(varargin{ilayer}{3}{2})
                        [mx mn] = slover('volmaxmin', so.img(ilayer).vol);
                        rng = [];
                        if isempty(varargin{ilayer}{3}{1})
                            rng(1) = mn;
                        else
                            rng(1) = varargin{ilayer}{3}{1};
                        end
                        if isempty(varargin{ilayer}{3}{2})
                            rng(2) = mx;
                        else
                            rng(2) = varargin{ilayer}{3}{2};
                        end
                        so.img(ilayer).range = rng;   
                    else
                        so.img(ilayer).range = [varargin{ilayer}{3}{1} varargin{ilayer}{3}{2}];
                    end
                case {'.mat'}
                     % Statspec
                     if isempty(varargin{ilayer}{3})                                   % If statspec is not specified
                        load(fl);                                                      % Load SPM.mat file in workspace
                        k = CorrClusTh(SPM,0.001,0.05);                                % Determine extent threshold for cluster-FWE correction
                        spc = {1,0.001,'none',k};

%                         [hReg,xSPM] = results_nogui('pv','spm',fl,'con',spc{1},'thr',spc{2},'ptype',spc{3},'ext',spc{4},'disp',0); %
                        [junk,xSPM] = results_nogui('pv','spm',fl,'con',spc{1},'thr',spc{2},'ptype',spc{3},'ext',spc{4},'disp',0);
                        if isempty(xSPM.Z) % If no suprathreshold voxels, use k = 0 (ie. P<.001 unc.)
                            clear hReg xSPM
%                             [hReg,xSPM] = results_nogui('pv','spm',fl,'con',spc{1},'thr',spc{2},'ptype',spc{3},'ext',0,'disp',0);
                            [junk,xSPM] = results_nogui('pv','spm',fl,'con',spc{1},'thr',spc{2},'ptype',spc{3},'ext',0,'disp',0);
                        end
                    else
                        if size(varargin{ilayer}{3},2) == 4
                            load(fl);
%                             [hReg,xSPM] = results_nogui('pv','spm',fl,'con',varargin{ilayer}{3}{1},'thr',varargin{ilayer}{3}{2},'ptype',varargin{ilayer}{3}{3},'ext',varargin{ilayer}{3}{4},'disp',0);
                            [junk,xSPM] = results_nogui('pv','spm',fl,'con',varargin{ilayer}{3}{1},'thr',varargin{ilayer}{3}{2},'ptype',varargin{ilayer}{3}{3},'ext',varargin{ilayer}{3}{4},'disp',0);
                        elseif size(varargin{ilayer}{3},2) == 8
                            load(fl);
%                             [hReg,xSPM] = results_nogui('pv','spm',fl,'con',varargin{ilayer}{3}{1},'thr',varargin{ilayer}{3}{2},'ptype',varargin{ilayer}{3}{3},'ext',varargin{ilayer}{3}{4},'mspm',varargin{ilayer}{3}{5},'mcon',varargin{ilayer}{3}{6},'mthr',varargin{ilayer}{3}{7},'mtype',varargin{ilayer}{3}{8},'disp',0);
                            [junk,xSPM] = results_nogui('pv','spm',fl,'con',varargin{ilayer}{3}{1},'thr',varargin{ilayer}{3}{2},'ptype',varargin{ilayer}{3}{3},'ext',varargin{ilayer}{3}{4},'mspm',varargin{ilayer}{3}{5},'mcon',varargin{ilayer}{3}{6},'mthr',varargin{ilayer}{3}{7},'mtype',varargin{ilayer}{3}{8},'disp',0);
                        elseif size(varargin{ilayer}{3},2) == 9
                            load(fl);
%                             [hReg,xSPM] = results_nogui('pv','spm',fl,'con',varargin{ilayer}{3}{1},'thr',varargin{ilayer}{3}{2},'ptype',varargin{ilayer}{3}{3},'ext',varargin{ilayer}{3}{4},'mspm',varargin{ilayer}{3}{5},'mcon',varargin{ilayer}{3}{6},'mthr',varargin{ilayer}{3}{7},'mtype',varargin{ilayer}{3}{8},'conj',varargin{ilayer}{3}{9},'disp',0);
                            [junk,xSPM] = results_nogui('pv','spm',fl,'con',varargin{ilayer}{3}{1},'thr',varargin{ilayer}{3}{2},'ptype',varargin{ilayer}{3}{3},'ext',varargin{ilayer}{3}{4},'mspm',varargin{ilayer}{3}{5},'mcon',varargin{ilayer}{3}{6},'mthr',varargin{ilayer}{3}{7},'mtype',varargin{ilayer}{3}{8},'conj',varargin{ilayer}{3}{9},'disp',0);
                     end
                    end
                     
                     % Header info of layer
                    so.img(ilayer).vol = slover('blobs2vol', xSPM.XYZ, xSPM.Z, xSPM.M);
                    so.img(ilayer).range = [xSPM.u max(xSPM.Z(:))];
                     
                otherwise
                    error('Wrong file format');
            end
            
            % Color map
            if isempty(varargin{ilayer}{4})
                so.img(ilayer).cmap = slover('getcmap','autumn'); % default
            else
                if ischar(varargin{ilayer}{4})
                    so.img(ilayer).cmap = slover('getcmap',varargin{ilayer}{4});
                elseif isnumeric(varargin{ilayer}{4})
                    so.img(ilayer).cmap = varargin{ilayer}{4};
                end
            end
                
            % Other parameters
            so.img(ilayer).outofrange = {0 64};
            so.img(ilayer).nancol = [0];
                
        case 'truecolor'
            
            % Layer type
            so.img(ilayer).type = 'truecolor';
            
            % Header info of layer
            if isempty(varargin{ilayer}{2})
                mesg = ['Select image file for layer ',num2str(ilayer)];
                so.img(ilayer).vol = spm_vol(spm_select(1,'image',mesg));
            else
%                 so.img(ilayer).vol = spm_vol(varargin{ilayer}{2});
                v_b = spm_vol(varargin{ilayer}{2});
                [y_b xyz_b] = spm_read_vols(v_b);
                so.img(ilayer).vol = slover('matrix2vol', y_b, v_b.mat);
            end

            % Range 
            if isempty(varargin{ilayer}{3})
                [mx mn] = slover('volmaxmin', so.img(ilayer).vol);
                so.img(ilayer).range = [mn mx];
            elseif isempty(varargin{ilayer}{3}{1}) | isempty(varargin{ilayer}{3}{2})
                [mx mn] = slover('volmaxmin', so.img(ilayer).vol);
                rng = [];
                if isempty(varargin{ilayer}{3}{1})
                    rng(1) = mn;
                else
                    rng(1) = varargin{ilayer}{3}{1};
                end
                if isempty(varargin{ilayer}{3}{2})
                    rng(2) = mx;
                else
                    rng(2) = varargin{ilayer}{3}{2};
                end
                so.img(ilayer).range = rng;                
            else
                so.img(ilayer).range = [varargin{ilayer}{3}{1} varargin{ilayer}{3}{2}];
            end

            % Color map
            if isempty(varargin{ilayer}{4})
                so.img(ilayer).cmap = slover('getcmap','actc'); % default
            else
                if ischar(varargin{ilayer}{4})
                    so.img(ilayer).cmap = slover('getcmap',varargin{ilayer}{4});
                elseif isnumeric(varargin{ilayer}{4})
                    so.img(ilayer).cmap = varargin{ilayer}{4};
                end
            end

            if size(varargin{ilayer},2) == 5
                so.img(ilayer).prop = varargin{ilayer}{5};
                y_b(find(y_b(:) == 1)) = NaN;                            % Replace ones with nans and load vol again
                so.img(ilayer).vol = slover('matrix2vol', y_b, v_b.mat); %  ... and load vol again
            else
%                 so.img(ilayer).prop = 1;
            end
            
            % Other parameters
            so.img(ilayer).outofrange = {0 64};
            so.img(ilayer).nancol = [0];
%             so.img(ilayer).outofrange = {0 64};
%             so.img(ilayer).nancol = [1 1 1];
            
        case 'contour'
            
            % Layer type
            so.img(ilayer).type = 'contour';
            
            % Header info of layer
            if isempty(varargin{1}{2})
                mesg = ['Select image file for layer ',num2str(ilayer)];
                so.img(ilayer).vol = spm_vol(spm_select(1,'image',mesg));
            else
                v_b = spm_vol(varargin{ilayer}{2});
                [y_b xyz_b] = spm_read_vols(v_b);
                so.img(ilayer).vol = slover('matrix2vol', y_b, v_b.mat);
            end

            % Line spec
            if isempty(varargin{ilayer}{3})
                so.img(ilayer).linespec = 'w-';
            else
                so.img(ilayer).linespec = varargin{ilayer}{3};
            end
            
            
            % Other parameters
            so.img(ilayer).outofrange = {0 64};
            so.img(ilayer).nancol = [0];
%             so.img(ilayer).outofrange = {0 64};
%             so.img(ilayer).nancol = [1 1 1];
            so.img(ilayer).linewidth = 0.5;
            

    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confgure general settings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orientation
if size(varargin{nargin},2) >= 1 & ~isempty(varargin{nargin}{1})
    so.transform = varargin{nargin}{1};
else
    so.transform = 'axial';
end

% Slices
if size(varargin{nargin},2) >= 2 & ~isempty(varargin{nargin}{2})
    so.slices = varargin{nargin}{2};
else
    so.slices = [-40 -20 0 20 40];
end

% Number of columns in figure
if size(varargin{nargin},2) >= 3 & ~isempty(varargin{nargin}{3})
    so.xslices = varargin{nargin}{3};
else
    so.xslices = 5;
end

% Color bars to display
if size(varargin{nargin},2) >= 4 & ~isempty(varargin{nargin}{4})
    so.cbar = varargin{nargin}{4};
else
    so.cbar = [];
end

% Labels to display
if size(varargin{nargin},2) >= 5 & ~isempty(varargin{nargin}{5})
    so.labels = varargin{nargin}{5};
else
    so.labels = struct('color','k');
end

% Resolution
if size(varargin{nargin},2) >= 6 & ~isempty(varargin{nargin}{6})
    so.printstr = ['print -dtiff -r',num2str(varargin{nargin}{6})];
else
    so.printstr = ['print -dtiff -r300'];
end

% Output filename
if size(varargin{nargin},2) >= 7 & ~isempty(varargin{nargin}{7})
    so.printfile = varargin{nargin}{7};
else
    so.printfile = fullfile(pwd,['fMRIfig_',datestr(now,'yyyymmmdd'),'.tiff']);
end

% Optional hard coded configurations
% so.area = struct('position',[0.025 0 0.9 1],'units','normalized','halign','left','valign','bottom');
so.area = struct('position',[0 0 1 1],'units','normalized','halign','left','valign','bottom');

% Fill empty fields
so = fill_defaults(so);

% Output filename
if size(varargin{nargin},2) >= 8 & ~isempty(varargin{nargin}{8})
    paperwidth = varargin{nargin}{8};
else
    paperwidth = 17;
end

if settings{9}
    
    % Configure figure
    X=1;Y=2;Z=3;
    dims = so.slicedef;
    xmm = dims(X,1):dims(X,2):dims(X,3);
    ymm = dims(Y,1):dims(Y,2):dims(Y,3);
    zmm = so.slices;
    [y x] = meshgrid(ymm,xmm');
    vdims = [length(xmm),length(ymm),length(zmm)];

    scrsz = get(0,'ScreenSize');
    sof = figure(so.figure);

    set(sof,'Units','centimeters');
    rows = ceil(vdims(3)/so.xslices);
    panelw = paperwidth/so.xslices;
    panelh = panelw*vdims(2)/vdims(1);
    figh =  rows*panelh;

    set(sof,'PaperUnits','centimeters');
    set(sof,'PaperSize',[paperwidth,figh]);
    set(sof,'PaperPosition',[0 0 paperwidth,figh]);
    set(sof,'Position',[0 scrsz(4)/2 paperwidth,figh]);
    set(sof,'Color',[1 1 1]);

    % Produce figure
    so = paint(so);
    % assignin('base','so',so);

    axh = findobj(sof,'XColor',[0 0 0]);        % Obtain handles of axes with XColor is black
    for iax = 1:size(axh,1)
        set(axh(iax),'XColor',[1 1 1],'YColor',[1 1 1]);  % Set axes to white
    end

    % axh = findobj(sof,'XColor',[1 1 1]);        % Obtain handles of axes with XColor is white
    % for iax = 1:size(axh,1)
    %     set(axh(iax),'XColor',[0 0 0],'YColor',[0 0 0]);  % Set axes to black
    % end

    % Print figure
    print_fig(so);
end
varargout{1} = so;