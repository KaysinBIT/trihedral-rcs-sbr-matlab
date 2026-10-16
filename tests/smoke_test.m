function smoke_test
% Quick reproduction and packaged-data consistency check.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'));
assert(strcmp(pattern_output_dir,fullfile(root,'results')));
a=43.6e-3; f=340e9; lambda=299792458/f;
r=sbr_pattern_exit(a,f,45,asind(1/sqrt(3)),384);
reference=4*pi*a^4/(3*lambda^2);
assert(abs(r.HH/reference-1)<1e-10);
assert(abs(r.VV/reference-1)<1e-10);
assert(all(isfinite([r.HH,r.VV,r.HV,r.VH])));
coarse=sbr_pattern_exit(a,f,26.5,24,768);
fine=sbr_pattern_exit(a,f,26.5,24,1536);
assert(abs(10*log10(coarse.HH/fine.HH))<0.03);
D=load(fullfile(root,'results','rcs_pattern_data.mat'));
assert(isequal(size(D.HH),[177,177]));
assert(all(isfinite(D.HH(:))) && all(D.HH(:)>0));
assert(all(isfinite(D.VV(:))) && all(D.VV(:)>0));
assert(D.summary.is_full_wave==false);
assert(strcmp(D.summary.sampling_method,...
    'Uniform exit-aperture sampling with reverse-ray visibility and coherent Huygens integration'));
T=readtable(fullfile(root,'results','rcs_azimuth_elevation.csv'));
assert(height(T)==177*177);
assert(max(abs(T.HH_m2-D.HH(:)))<1e-10);
assert(max(abs(T.VV_m2-D.VV(:)))<1e-10);
fprintf('PASS: portable paths, PEC reference, HH/VV symmetry, grazing-branch stability, MAT/CSV consistency.\n');
fprintf('Boresight HH %.9f dBsm; packaged angular points %d.\n',10*log10(r.HH),height(T));
end
