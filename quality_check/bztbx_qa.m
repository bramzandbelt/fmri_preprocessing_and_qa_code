function bztbx_qa(varargin)
% _________________________________________________________________________
% QUALITY ASSURANCE Performs a basic quality check of fMRI time series
%
% FORMAT bztbx_qa('preproc',files,t1img,rpfile,isi,hpf,outputdir)
% DESCRIPTION     For analysis of realigned scans
% files           - filenames of realigned scans (char array)
% t1img           - filename of coregistered T1-weighted scan (char array)
% rpfile          - filename of realignment params (char array)
% isi             - interscan interval (scalar, in seconds)
% hpf             - highpass filter cutoff (scalar, in seconds)
% outputdir       - directory where files should be written (char array)
% 
% ________________________________________________________________________
% Bram B. Zandbelt, July 2008
% b.b.zandbelt@umcutrecht.nl
%
% version 1.0

% Defining inputs
if nargin == 0
    error('No input arguments provided. Type help qualitycheck for usage')
end

switch lower(varargin{1})
    case 'preproc'
        if nargin  ~= 7
            error('Wrong number of input arguments. Type help qualitycheck for usage')
        end
        ts = varargin{2};
        t1img = varargin{3};
        rpfile = varargin{4};
        isi = varargin{5};
        hp_cutoff = varargin{6};
        outputdir = varargin{7};
    case '1stlevstat'
        if nargin  ~= 2
            error('Wrong number of input arguments. Type help qualitycheck for usage')
        end
        spmmat = varargin{2};
end

switch lower(varargin{1})
    case 'preproc'

        if isempty(ts) == 1
            ts = spm_select(Inf,'image','Select realigned fMRI data');
        end

        if isempty(t1img) == 1
            t1img = spm_select(1,'image','Select T1-weighted scan');
        end


        if isempty(rpfile) == 1
            rpfile = spm_select(1,'any','Select realignment parameters file',[],pwd,'^rp.*\.txt$');
        end

        if isempty(outputdir) == 1
            outputdir = spm_select(1,'dir','Select directory where output should be written');
        end

        if exist(outputdir,'dir') ~= 7
            mkdir (outputdir)
        end

        preswd = pwd; % Present working directory
        
        [p n] = spm_fileparts(deblank(ts(1,:)));
        
        
        
        % Calculating mean raw ts, mean hpf ts, tmean, tsd, tsnr
        [maskmatrix,mask] = bztbx_automask('robust',t1img,deblank(ts(1,:)));
        clear maskmatrix
        
        bztbx_tsnr('hpf',ts,isi,hp_cutoff,outputdir,mask);

        % Loading raw ts data
        tsraw = spm_select('FPList',outputdir,['^tsraw_',n,'.*\.mat$']);
        load(tsraw);
        
        % Loading hpf ts data
        tshpf = spm_select('FPList',outputdir,['^tshpf_',n,'.*\.mat$']);
        load(tshpf);
        
        % Loading global psc ts data
        tspsc = spm_select('FPList',outputdir,['^tspsc_',n,'.*\.mat$']);
        load(tspsc);
        
        % Loading local psc ts data
        tspscsl = spm_select('FPList',outputdir,['^tspscsl_',n,'.*\.mat$']);
        load(tspscsl);
        
        % Loading realignment parameter file
        rpmat = load(rpfile);
        
%         Generating figure
        fig = spm_figure('FindWin','Graphics');
        if isempty(fig)
            fig = spm_figure('Create','Graphics');
        end;
        set(0,'CurrentFigure',fig);
        spm_figure('Clear','Graphics');
% 
%         Calculate the number of bins required
        nScan = size(ts,1);
        bins = 30;                                                  % Approximate maximum number of X ticks
        minbin = 10;                                                % Minimum distance between X ticks
        factor = ceil(nScan/(bins*minbin));                        % Factor to calculate X tick vector
        xtickvec = factor*minbin:factor*minbin:nScan;            % X tick vector

        %Plot 1: raw and high-pass filtered time series
        subplot(8,1,1)                                              
        ax = plotyy(1:nScan,tsraw,1:nScan,tshpf);
        
        title({[datestr(now,'dd-mmm-yyyy HH:MM PM'),': Quality check for ',n(1:end-5)];},...
                    'FontWeight','bold',...
                    'Interpreter','none'...
                    );
        grid on
        
        set(get(ax(1),'Ylabel'),'String','S (a.u.)');
        set(get(ax(2),'Ylabel'),'String','S (a.u.)');

        set(ax(1),'XTick',xtickvec);
        set(ax(2),'XTick',xtickvec);

        set(ax(1),'XLim',[0,nScan]);
        set(ax(2),'XLim',[0,nScan]);
        
        figtitle{1} = ['raw & high-pass (cutoff: ',num2str(hp_cutoff),' s) filtered signal'];
        axs = get(ax);
        x1 = 0.02*(axs(1).XLim(2)-axs(1).XLim(1))+axs(1).XLim(1);
        y1 = 0.75*(axs(1).YLim(2)-axs(1).YLim(1))+axs(1).YLim(1);
        text(x1,y1,figtitle{1},'Color','r','Interpreter','none')
        
        %Plot 2: global percent signal change
        subplot(8,1,2)
        plot(tspsc)
        
        
        
        grid on
        axis tight
        ylabel('dS (%)');
        hgca(2) = get(gca);
        figtitle{2} = ['global % change in HPF (',num2str(hp_cutoff),' s) signal'];
        x1 = 0.02*(hgca(2).XLim(2)-hgca(2).XLim(1))+hgca(2).XLim(1);
        y1 = 0.75*(hgca(2).YLim(2)-hgca(2).YLim(1))+hgca(2).YLim(1);
        text(x1,y1,figtitle{2},'Color','r','Interpreter','none')
        set(gca,'XTick',xtickvec);

        % Plot 3: local (i.e. per slice) percent signal change
        subplot(8,1,3:4)                                          
        spm_figure('ColorMap','jet')

        clim = [0 2]; % Color bar limits
        
        imagesc(tspscsl',clim);
        colorbar('location','NorthOutside')

        grid on
        axis tight
        ylabel('slice in dim 3');
        hgca(3) = get(gca);
        figtitle{3} = ['local % change in HPF (',num2str(hp_cutoff),' s) signal per slice'];
        x1 = 0.02*(hgca(3).XLim(2)-hgca(3).XLim(1))+hgca(3).XLim(1);
        y1 = 0.25*(hgca(3).YLim(2)-hgca(3).YLim(1))+hgca(3).YLim(1);
        text(x1,y1,figtitle{3},'Color','w','Interpreter','none')
        set(gca,'XTick',xtickvec);
        set(gca,'XLim',[0,nScan]);

        % Plot 4: realignment parameters
        [rpp rpn] = spm_fileparts(deblank(ts(1,:)));
        
        subplot(8,1,5:6)    

        [ax,h1,h2] = plotyy(1:nScan,rpmat(:,1:3),1:nScan,rpmat(:,4:6));
        grid on

        set(h1,'LineStyle','-')
        set(h2,'LineStyle','--')

        set(h1(1),'Color',[0 0 1]);
        set(h1(2),'Color',[0 0.5 0]);
        set(h1(3),'Color',[1 0 0]);

        set(h2(1),'Color',[0 0 1]);
        set(h2(2),'Color',[0 0.5 0]);
        set(h2(3),'Color',[1 0 0]);

        set(get(ax(1),'Ylabel'),'String','Translation (mm)');
        set(get(ax(2),'Ylabel'),'String','Rotation (deg)');

        set(ax(1),'XTick',xtickvec);
        set(ax(2),'XTick',xtickvec);

        set(ax(1),'XLim',[0,nScan]);
        set(ax(2),'XLim',[0,nScan]);

        axs = get(ax);
        figtitle{4} = [rpn,' - translations & rotations'];
        x1 = 0.02*(axs(1).XLim(2)-axs(1).XLim(1))+axs(1).XLim(1);
        y1 = 0.75*(axs(1).YLim(2)-axs(1).YLim(1))+axs(1).YLim(1);
        text(x1,y1,figtitle{4},'Color','r','Interpreter','none');

        legend('location','NorthOutside','orientation','horizontal');
        [legh,objh,outh] = legend;
        set(outh(1),'DisplayName','x & pitch');
        set(outh(2),'DisplayName','y & roll');
        set(outh(3),'DisplayName','z & yaw');
        set(legh,'Box','off')

%         Plot of 'fast motion'
        subplot(8,1,7:8) 

        rpmatdiff = rpmat(2:end,:) - rpmat(1:end-1,:);
        r = 65; % radius (i.e. distance between cortex and center of brain, in mm)
        x = [-r 0 0 1]';
        fastmotion = zeros(size(rpmatdiff,1),1);

        for idyn = 1:size(rpmatdiff,1)
            q = rpmatdiff(idyn,:);
            q(4:6) = q(4:6)*pi/180;     % convert degrees into radians
            T = [1 0 0 q(1);0 1 0 q(2);0 0 1 q(3); 0 0 0 1];    % setting up the translation matrix
            % setting up the rotation matrices
            R1 = [1 0 0 0;0 cos(q(4)) sin(q(4)) 0; 0 -sin(q(4)) cos(q(4)) 0; 0 0 0 1];
            R2 = [cos(q(5)) 0 sin(q(5)) 0; 0 1 0 0; -sin(q(5)) 0 cos(q(5)) 0; 0 0 0 1];
            R3 = [cos(q(6)) sin(q(6)) 0 0;-sin(q(6)) cos(q(6)) 0 0; 0 0 1 0;0 0 0 1];
            R = R1*R2*R3;
            M = T*R;
            y = M*x;
            fastmotion(idyn) = sqrt((y(1)-x(1))^2 + (y(2)-x(2))^2 + (y(3)-x(3))^2);

            clear q T R* M y
        end

        plot(2:nScan,fastmotion);
        grid on
        xlabel('scan')
        ylabel('motion (mm)');

        hgca(5) = get(gca);
        figtitle{5} = [rpn,' - motion per scan (estimated @ ',num2str(r),' mm from brain center) '];
        x1 = 0.02*(hgca(5).XLim(2)-hgca(5).XLim(1))+hgca(5).XLim(1);
        y1 = 0.75*(hgca(5).YLim(2)-hgca(5).YLim(1))+hgca(5).YLim(1);
        text(x1,y1,figtitle{5},'Color','r','Interpreter','none')
        set(gca,'XTick',xtickvec);
        set(gca,'XLim',[0,nScan]);


        % Save percent signal change (PSC) and fast movement (FMV) data
%         PSC.mean        = mean(tspsc);
%         PSC.std         = std(tspsc);
%         PSC.min         = min(tspsc);
%         PSC.max         = max(tspsc);
%         PSC.CI95_low    = PSC.mean - tinv(.975,size(tspsc,1) - 1)*std(tspsc)/sqrt(size(tspsc,1));
%         PSC.CI95_high   = PSC.mean + tinv(.975,size(tspsc,1) - 1)*std(tspsc)/sqrt(size(tspsc,1));

        FMV.mean        = mean(fastmotion);
        FMV.std         = std(fastmotion);
        FMV.min         = min(fastmotion);
        FMV.max         = max(fastmotion);
        FMV.CI95_low    = FMV.mean - tinv(.975,size(fastmotion,1) - 1)*nanstd(fastmotion)/sqrt(size(fastmotion,2));
        FMV.CI95_high   = FMV.mean + tinv(.975,size(fastmotion,1) - 1)*nanstd(fastmotion)/sqrt(size(fastmotion,2));

%         txtrep = horzcat(PSC.mean,PSC.std,PSC.min,PSC.max,PSC.CI95_low,PSC.CI95_high, ...
%                          FMV.mean,FMV.std,FMV.min,FMV.max,FMV.CI95_low,FMV.CI95_high);

        txtrep = horzcat(FMV.mean,FMV.std,FMV.min,FMV.max,FMV.CI95_low,FMV.CI95_high);
             
                     
        fid = fopen(fullfile(outputdir,['motiondata_',n(1:end-5),'.txt']),'w+');
        fprintf(fid,'FMV_Mean \t FMV_Std \t FMV_Min \t FMV_Max \t FMV_95CI_low \t FMV_95CI_high \n');
        fprintf(fid,'%0.4f    \t %0.4f   \t %0.4f   \t %0.4f   \t %0.4f        \t %0.4f', txtrep);
        fprintf(fid,'\n');

        % Save figure
        spm_print;
        spm_figure('ColorMap','gray');


        % Temporal Mean, SD and SNR data
        
        for imap = 1:3

            fig = spm_figure('FindWin','Graphics');
            if isempty(fig)
                fig = spm_figure('Create','Graphics');
            end;
            set(0,'CurrentFigure',fig);
            spm_figure('Clear','Graphics');

            SO = slover();

            if imap == 1
                overlay = spm_select('FPList',outputdir,['^tmean_',n,'.*\.nii$']);
            elseif imap == 2
                overlay = spm_select('FPList',outputdir,['^tsd_',n,'.*\.nii$']);
            elseif imap == 3
                overlay = spm_select('FPList',outputdir,['^tsnr_',n,'.*\.nii$']);
            end

            % Slice overlay settings
            SO.img(1).vol = spm_vol(t1img);                 % T1-weighted scan
            SO.img(2).vol = spm_vol(overlay);               % Overlay

            if imap == 1
                SO.img(2).range = [0 250];                 % Color bar range of overlay
            elseif imap == 2
                SO.img(2).range = [0 10];                 % BBZ: Sep 8, 2012 - Adjust for Vanderbilt macaque data
%                 SO.img(2).range = [0 2];
            elseif imap == 3
                SO.img(2).range = [0 100];               % BBZ: Sep 8, 2012 - Adjust for Vanderbilt macaque data
%                 SO.img(2).range = [0 250];
            end

            [cmap] = slover('getcmap', 'actc');     % Create color map
            SO.img(2).cmap = cmap;                          
            SO.cbar = 2;                                    % Color bar represents values in map 2
            SO.transform = 'axial';                         % Slice orientation
            SO.slices = 0:4:40;                           % BBZ: Sep 8, 2012 - Adjust for macaque data
%             SO.slices = -16:8:64;                           % Coordinates of slices
            SO.printstr = 'print -noui -painters ';         % Print settings
            paint(SO);                                 % Make overlay

            % Editing figure header
            hs = getappdata(gcf,'PrintHeaderHeaderSpec');
            if isempty(hs)
             hs = struct('dateformat','none',...
             'string','',...
             'fontname','Arial',...
             'fontsize',12,... % in points
             'fontweight','normal',...
             'fontangle','normal',...
             'margin',72); % in points
            end

            %currentPage = 1;
            %totalPages = 10;
            if imap == 1
                hs.string = [datestr(now,'dd-mmm-yyyy HH:MM PM'),': tmean map from tsnr.m for ',n(1:end-5)];
            elseif imap == 2
                hs.string = [datestr(now,'dd-mmm-yyyy HH:MM PM'),': tsd map from tsnr.m for ',n(1:end-5)];
            elseif imap == 3
                hs.string = [datestr(now,'dd-mmm-yyyy HH:MM PM'),': tsnr map from tsnr.m for ',n(1:end-5)];
            end
            setappdata(gcf,'PrintHeaderHeaderSpec',hs);

            set(gcf,'PaperPositionMode','auto');


            spm_print;

%             % Added 25-Feb-2011
%             SO.transform = 'sagittal';                         % Slice orientation
%             SO.slices = -64:16:64;    
%             paint(SO);
%             spm_print
%             
%             clear overlay SO cmap warnstr

            
            
            
            
        end
   
    case '1stlevstat'
        for ifls = 1:size(spmmat,1)
            [p] = spm_fileparts(strtrim(spmmat(ifls,:)));
            cd(p)
            resms = spm_select('FPList',p,'^ResMS.*img$');
            rpv = spm_select('FPList',p,'^RPV.*img$');
            out1 = 'ResStd.img';
            out2 = 'FWHM.img';
            eq1 = 'sqrt(i1)';
            eq2 = 'i1.^(-1/3)';
            flags = {0 0 4 1};
            spm_imcalc_ui(resms,out1,eq1,flags);
            spm_imcalc_ui(rpv,out2,eq2,flags);
            
            clear p n e v resms rpv
        end
        
        
end
