function run_pattern_sweep
% 340 GHz azimuth/elevation pattern using the previously verified PEC SBR.
src=fileparts(mfilename('fullpath')); addpath(src);
outdir=pattern_output_dir; if ~isfolder(outdir), mkdir(outdir); end
diary(fullfile(outdir,'sweep.log')); diaryCleanup=onCleanup(@()diary('off'));
a=43.6e-3; f=340e9; el0=asind(1/sqrt(3));
az=1:0.5:89; el=1:0.5:89; Nmap=384; Ncut=768;
[AZ,EL]=meshgrid(az,el); aa=AZ(:); ee=EL(:); total=numel(aa);
values=nan(total,5); checkpoint=fullfile(outdir,'map_exit_checkpoint.mat');
if isfile(checkpoint)
    old=load(checkpoint);
    assert(isequal(old.az,az) && isequal(old.el,el) && old.Nmap==Nmap && old.a==a && old.f==f);
    values=old.values;
end
pool=gcp('nocreate'); if isempty(pool), pool=parpool('local',4); end
poolCleanup=onCleanup(@()delete(pool));
fprintf('MATLAB %s; PEC a=%.6g m; f=%.6g Hz\n',version,a,f);
fprintf('Map %d x %d = %d dual-pol angles; N=%d (%d rays/angle)\n',...
    numel(el),numel(az),total,Nmap,Nmap^2);
clock0=tic;
missing=find(isnan(values(:,1)));
for first=1:1500:numel(missing)
    ids=missing(first:min(first+1499,numel(missing)));
    block=zeros(numel(ids),5);
    parfor j=1:numel(ids)
        idx=ids(j);
        r=sbr_pattern_exit(a,f,aa(idx),ee(idx),Nmap);
        block(j,:)=[r.HH,r.VV,r.HV,r.VH,r.tripleHH];
    end
    values(ids,:)=block;
    save(checkpoint,'values','az','el','Nmap','a','f');
    fprintf('Map completed %d/%d; elapsed %.1f s\n',nnz(~isnan(values(:,1))),total,toc(clock0));
end
mapSeconds=toc(clock0);
HH=reshape(values(:,1),size(AZ)); VV=reshape(values(:,2),size(AZ));
HV=reshape(values(:,3),size(AZ)); VH=reshape(values(:,4),size(AZ));
tripleHH=reshape(values(:,5),size(AZ));
boresight=sbr_pattern_exit(a,f,45,el0,1536);

cutAz=0.25:0.25:89.75;
cutEl=sort(unique([0.25:0.25:89.75,el0]));
cutAzValues=zeros(numel(cutAz),2); cutElValues=zeros(numel(cutEl),2);
parfor i=1:numel(cutAz)
    r=sbr_pattern_exit(a,f,cutAz(i),el0,Ncut);
    cutAzValues(i,:)=[r.HH,r.VV];
end
parfor i=1:numel(cutEl)
    r=sbr_pattern_exit(a,f,45,cutEl(i),Ncut);
    cutElValues(i,:)=[r.HH,r.VV];
end
fprintf('Principal cuts completed: %d azimuth + %d elevation points.\n',numel(cutAz),numel(cutEl));
save(fullfile(outdir,'cuts_checkpoint.mat'),'cutAz','cutEl','cutAzValues','cutElValues','Ncut');

% Compare map ray density against 4x and 16x as many rays at 90 angles.
[VA,VE]=meshgrid([1,5,15,30,45,60,75,85,89],[1,5,15,25,el0,45,60,75,85,89]);
vAz=[VA(:);26.5;63.5;45;45;45;44.5;45.5];
vEl=[VE(:);24;24;56.5;35;35.5;35.5;35.5]; validation=zeros(numel(vAz),6);
validationFile=fullfile(outdir,'ray_density_validation_exit.csv');
if isfile(validationFile)
    cached=readtable(validationFile);
    assert(max(abs(cached.Azimuth_deg-vAz))<1e-10 && max(abs(cached.Elevation_deg-vEl))<1e-10);
    validation=cached{:,3:8};
else
    parfor i=1:numel(vAz)
        r1=sbr_pattern_exit(a,f,vAz(i),vEl(i),384);
        r2=sbr_pattern_exit(a,f,vAz(i),vEl(i),768);
        r3=sbr_pattern_exit(a,f,vAz(i),vEl(i),1536);
        validation(i,:)=10*log10([r1.HH,r2.HH,r3.HH,r1.VV,r2.VV,r3.VV]);
    end
end
validationTable=array2table([vAz,vEl,validation],...
    'VariableNames',{'Azimuth_deg','Elevation_deg','HH_N384_dBsm','HH_N768_dBsm',...
    'HH_N1536_dBsm','VV_N384_dBsm','VV_N768_dBsm','VV_N1536_dBsm'});
validationTable.HH_MapMinusFine_dB=validation(:,1)-validation(:,3);
validationTable.VV_MapMinusFine_dB=validation(:,4)-validation(:,6);
writetable(validationTable,validationFile);
referencePeak=10*log10(boresight.HH);
strongHH=validation(:,3)>=referencePeak-10;
strongVV=validation(:,6)>=referencePeak-10;
mainHH=validation(:,3)>=referencePeak-3;
mainVV=validation(:,6)>=referencePeak-3;
moderateHH=validation(:,3)>=referencePeak-20;
moderateVV=validation(:,6)>=referencePeak-20;
validationStats=struct(...
    'tested_angles',numel(vAz),...
    'max_map_error_above_peak_minus_3dB',max(abs([validation(mainHH,1)-validation(mainHH,3);validation(mainVV,4)-validation(mainVV,6)])),...
    'max_map_error_above_peak_minus_10dB',max(abs([validation(strongHH,1)-validation(strongHH,3);validation(strongVV,4)-validation(strongVV,6)])),...
    'max_map_error_above_peak_minus_20dB',max(abs([validation(moderateHH,1)-validation(moderateHH,3);validation(moderateVV,4)-validation(moderateVV,6)])),...
    'max_tested_map_error_dB',max(abs([validation(:,1)-validation(:,3);validation(:,4)-validation(:,6)])));
fprintf('Ray density validation: strong max %.6f dB; moderate max %.6f dB; all max %.6f dB\n',...
    validationStats.max_map_error_above_peak_minus_10dB,...
    validationStats.max_map_error_above_peak_minus_20dB,validationStats.max_tested_map_error_dB);

[gridMax,index]=max(HH(:));
azWidth=half_power_width(cutAz,10*log10(cutAzValues(:,1)),45,referencePeak-3);
elWidth=half_power_width(cutEl,10*log10(cutElValues(:,1)),el0,referencePeak-3);
summary=struct('status','completed_PEC_SBR_high_frequency_approximation',...
    'sampling_method','Uniform exit-aperture sampling with reverse-ray visibility and coherent Huygens integration',...
    'frequency_Hz',f,'inner_leg_m',a,'matlab_version',version,...
    'azimuth_limits_deg',[az(1),az(end)],'elevation_limits_deg',[el(1),el(end)],...
    'map_angle_step_deg',0.5,'cut_angle_step_deg',0.25,...
    'map_aperture_samples_per_angle',Nmap^2,'cut_aperture_samples_per_angle',Ncut^2,...
    'max_exit_branches_per_aperture_sample',4,...
    'map_dual_polarization_angles',total,'map_elapsed_seconds_this_run',mapSeconds,...
    'boresight_az_deg',45,'boresight_el_deg',el0,...
    'boresight_HH_dBsm',10*log10(boresight.HH),'boresight_VV_dBsm',10*log10(boresight.VV),...
    'grid_peak_HH_dBsm',10*log10(gridMax),'grid_peak_az_deg',AZ(index),'grid_peak_el_deg',EL(index),...
    'HH_azimuth_3dB',azWidth,'HH_elevation_3dB',elWidth,...
    'ray_density_validation',validationStats,...
    'angular_domain','Open first octant: s_x,s_y,s_z>0; exterior/back directions are not calculated',...
    'display_floor_dBsm',-25,'is_full_wave',false,...
    'limitations',{{'Ideal zero-thickness PEC plates with exact orthogonality',...
    'No explicit UTD/PTD physical-edge diffraction, loss, roughness or assembly errors',...
    'Deep minima and near-grazing details are more sensitive to ray density',...
    'Ray-grid agreement does not bound error relative to a full-wave solution or measurement'}});
save(fullfile(outdir,'rcs_pattern_data.mat'),'a','f','az','el','HH','VV','HV','VH','tripleHH',...
    'cutAz','cutEl','cutAzValues','cutElValues','boresight','summary','validationTable');
T=array2table([aa,ee,10*log10(values(:,1:2)),values(:,1:2)],...
    'VariableNames',{'Azimuth_deg','Elevation_deg','HH_dBsm','VV_dBsm','HH_m2','VV_m2'});
writetable(T,fullfile(outdir,'rcs_azimuth_elevation.csv'));
writetable(array2table([cutAz',cutAz'-45,10*log10(cutAzValues)],...
    'VariableNames',{'Azimuth_deg','OffsetFromBoresight_deg','HH_dBsm','VV_dBsm'}),...
    fullfile(outdir,'azimuth_cut.csv'));
writetable(array2table([cutEl',cutEl'-el0,10*log10(cutElValues)],...
    'VariableNames',{'Elevation_deg','OffsetFromBoresight_deg','HH_dBsm','VV_dBsm'}),...
    fullfile(outdir,'elevation_cut.csv'));
fid=fopen(fullfile(outdir,'summary.json'),'w','n','UTF-8');
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true)); fclose(fid);
plot_pattern_sweep(outdir);
run_boresight_cuts;
fprintf('DONE: grid peak %.6f dBsm at az %.3f el %.3f deg; HH widths %.3f / %.3f deg.\n',...
    summary.grid_peak_HH_dBsm,summary.grid_peak_az_deg,summary.grid_peak_el_deg,...
    azWidth.width_deg,elWidth.width_deg);
end

function info=half_power_width(x,y,center,level)
x=x(:); y=y(:);
left=find(x<center & y<=level,1,'last');
right=find(x>center & y<=level,1,'first');
assert(~isempty(left) && ~isempty(right));
xl=x(left)+(level-y(left))*(x(left+1)-x(left))/(y(left+1)-y(left));
xr=x(right-1)+(level-y(right-1))*(x(right)-x(right-1))/(y(right)-y(right-1));
info=struct('left_deg',xl,'right_deg',xr,'width_deg',xr-xl,'threshold_dBsm',level);
end
