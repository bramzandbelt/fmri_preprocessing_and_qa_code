function bztbx_multireg(varargin)
% Multiple regression file construction
% FORMAT bztbx_multireg('default',srcfls,filename)
%
% FORMAT bztbx_multireg('ppi',srcfls,filename)
%

type = varargin{1};
files = varargin{2};
name = varargin{3};

[p n e v] = spm_fileparts(strtrim(files(1,:)));

switch type
    case 'default' % i.e. RSfMRI
        
        regressor = [];
        for ifile = 1:size(files,1);

            [regressor{ifile}.p,regressor{ifile}.n,regressor{ifile}.e,regressor{ifile}.v] = spm_fileparts(strtrim(files(ifile,:)));
            regressor{ifile}.signal = importdata(strtrim(files(ifile,:)));
        end
        
        multireg = [];
        multireg = regressor{1}.signal;
        if size(files,1) > 1

            for ifile = 2:size(files,1);
                multireg = horzcat(multireg,regressor{ifile}.signal);
            end

        end

        save(fullfile(regressor{1}.p,[name,'.txt']),'multireg','-ASCII');

    case 'ppi'
        
        
        if size(files,1) == 1
            
            PPI1 = load(files);
            if isempty(name)
                multiregstr = ['Multireg_PPI_' PPI1.PPI.name];
            else
                multiregstr = name;
            end
            
            multireg = horzcat(PPI1.PPI.ppi,PPI1.PPI.P,PPI1.PPI.Y);

            if spm_matlab_version_chk('7') >= 0,
                save(fullfile(p,[multiregstr,'.txt']),'-V6','-ASCII','multireg')
            else
                save(fullfile(p,[multiregstr,'.txt']),'-ASCII','multireg')
            end
            
        elseif size(files,1) == 2
            
            PPI1 = load(strtrim(files(1,:)));
            PPI2 = load(strtrim(files(2,:)));
            
            if isempty(name)
                multiregstr = ['Multireg_PPI_' PPI1.PPI.name '_and_' PPI2.PPI.name];
            else
                multiregstr = name;
            end
            
            multireg = horzcat(PPI1.PPI.ppi,PPI2.PPI.ppi,[PPI1.PPI.P-PPI2.PPI.P],PPI1.PPI.Y);
            
            if spm_matlab_version_chk('7') >= 0,
                save(fullfile(p,[multiregstr,'.txt']),'-V6','-ASCII','multireg')
            else
                save(fullfile(p,[multiregstr,'.txt']),'-ASCII','multireg')
            end
        
        end
end
