function verify_pattern_peak
% Recheck all highest sampled HH/VV candidates before accepting the plot.
src=fileparts(mfilename('fullpath')); outdir=pattern_output_dir;
file=fullfile(outdir,'rcs_pattern_data.mat'); D=load(file);
[AZ,EL]=meshgrid(D.az,D.el);
[~,ih]=sort(D.HH(:),'descend'); [~,iv]=sort(D.VV(:),'descend');
ids=unique([ih(1:12);iv(1:12)]); rows=zeros(numel(ids),6);
for i=1:numel(ids)
    id=ids(i); fine=sbr_pattern_exit(D.a,D.f,AZ(id),EL(id),1536);
    rows(i,:)=[AZ(id),EL(id),10*log10(D.HH(id)),10*log10(D.VV(id)),10*log10(fine.HH),10*log10(fine.VV)];
end
T=array2table(rows,'VariableNames',{'Azimuth_deg','Elevation_deg','HH_Map_dBsm','VV_Map_dBsm','HH_Fine_dBsm','VV_Fine_dBsm'});
writetable(T,fullfile(outdir,'peak_validation.csv'));
error=max(abs(rows(:,3:4)-rows(:,5:6)),[],'all');
assert(error<0.15,'Highest map candidates did not pass ray-density validation.');
D.summary.peak_validation=struct('candidate_angles',numel(ids),...
    'fine_aperture_samples_per_angle',1536^2,'max_change_dB',error,...
    'max_fine_HH_dBsm',max(rows(:,5)),'max_fine_VV_dBsm',max(rows(:,6)));
% Normalize metadata from the initial running version of this script.
if isfield(D.summary,'map_rays_per_angle')
    D.summary.map_aperture_samples_per_angle=D.summary.map_rays_per_angle;
    D.summary=rmfield(D.summary,'map_rays_per_angle');
end
if isfield(D.summary,'cut_rays_per_angle')
    D.summary.cut_aperture_samples_per_angle=D.summary.cut_rays_per_angle;
    D.summary=rmfield(D.summary,'cut_rays_per_angle');
end
D.summary.max_exit_branches_per_aperture_sample=4;
save(file,'-struct','D');
fid=fopen(fullfile(outdir,'summary.json'),'w','n','UTF-8');
fprintf(fid,'%s\n',jsonencode(D.summary,'PrettyPrint',true)); fclose(fid);
fprintf('Highest-peak recheck passed: %d angles, maximum change %.6g dB.\n',numel(ids),error);
end
