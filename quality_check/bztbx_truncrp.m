function bztbx_truncrp(varargin)
%
%
%

rpfile = varargin{1};
remove = varargin{2};       % how many time points to remove

[p n e v] = spm_fileparts(rpfile);
rp = load(rpfile);
rp_trunc = rp([1+remove:end],:);
save(fullfile(p,['t',n,e,v]),'rp_trunc','-ASCII');