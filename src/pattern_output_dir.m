function outdir=pattern_output_dir
% Portable location for published and regenerated simulation artifacts.
root=fileparts(fileparts(mfilename('fullpath')));
outdir=fullfile(root,'results');
if ~isfolder(outdir), mkdir(outdir); end
end
