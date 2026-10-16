function run_boresight_cuts
% Physical angular offsets about the aperture normal, not merely phi-45 deg.
src=fileparts(mfilename('fullpath')); outdir=pattern_output_dir;
addpath(src); D=load(fullfile(outdir,'rcs_pattern_data.mat'));
forward=[1,1,1]/sqrt(3); horizontal=[-1,1,0]/sqrt(2);
vertical=[-1,-1,2]/sqrt(6);
alpha=(-39:0.25:39)'; beta=D.cutEl(:)-D.summary.boresight_el_deg;
localHHVV=zeros(numel(alpha),2); globalAngles=zeros(numel(alpha),2);
for i=1:numel(alpha)
    s=cosd(alpha(i))*forward+sind(alpha(i))*horizontal;
    az=atan2d(s(2),s(1)); el=asind(s(3));
    r=sbr_pattern_exit(D.a,D.f,az,el,768);
    localHHVV(i,:)=[r.HH,r.VV]; globalAngles(i,:)=[az,el];
end
azdB=10*log10(localHHVV); eldB=10*log10(D.cutElValues);
threshold=D.summary.boresight_HH_dBsm-3;
azWidth=find_width(alpha,azdB(:,1),threshold);
elWidth=find_width(beta,eldB(:,1),threshold);
localSummary=struct('forward_axis',forward,'horizontal_axis',horizontal,'vertical_axis',vertical,...
    'direction_definition','s=cos(beta)*cos(alpha)*forward + cos(beta)*sin(alpha)*horizontal + sin(beta)*vertical',...
    'HH_azimuth_3dB',azWidth,'HH_elevation_3dB',elWidth,...
    'azimuth_cut_beta_deg',0,'elevation_cut_alpha_deg',0,'aperture_samples_per_angle',768^2);
save(fullfile(outdir,'boresight_local_cuts.mat'),'alpha','beta','localHHVV','globalAngles','localSummary','eldB');
T=array2table([alpha,globalAngles,azdB],...
    'VariableNames',{'LocalAzimuth_deg','GlobalAzimuth_deg','GlobalElevation_deg','HH_dBsm','VV_dBsm'});
writetable(T,fullfile(outdir,'boresight_azimuth_cut.csv'));
fid=fopen(fullfile(outdir,'boresight_cut_summary.json'),'w','n','UTF-8');
fprintf(fid,'%s\n',jsonencode(localSummary,'PrettyPrint',true)); fclose(fid);
set(groot,'defaultFigureVisible','off','defaultAxesFontName','Microsoft YaHei',...
    'defaultTextFontName','Microsoft YaHei','defaultAxesFontSize',12);
fig=figure('Color','w','Position',[50,50,1250,570]);
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
blue=[0.05,0.35,0.7]; orange=[0.82,0.32,0.08];
for plane=1:2
    nexttile;
    if plane==1
        x=alpha; y=azdB; w=azWidth; name='方位';
    else
        x=beta; y=eldB; w=elWidth; name='俯仰';
    end
    p1=plot(x,y(:,1),'-','LineWidth',1.6,'Color',blue); hold on;
    p2=plot(x,y(:,2),'--','LineWidth',1.4,'Color',orange);
    yline(threshold,':','−3 dB','LabelHorizontalAlignment','left','LineWidth',1.2);
    xline(0,':','LineWidth',1); grid on; box on;
    plot([w.left_deg,w.right_deg],[threshold,threshold],'ko','MarkerFaceColor','w','MarkerSize',5);
    xlabel([name '角：正对开口为 0°']); ylabel('单站 RCS (dBsm)');
    title(sprintf('%s切面：HH −3 dB 宽度 %.2f°',name,w.width_deg));
    xlim([min(x),max(x)]); ylim([-35,15]);
    legend([p1,p2],{'HH','VV'},'Location','south');
end
title(tl,{'340 GHz 三面角 RCS：以开口法向为中心的角度切面',...
    '内棱长 43.6 mm，PEC SBR；步长 0.25°；范围受第一卦限限制'},...
    'FontSize',14,'FontWeight','normal');
exportgraphics(fig,fullfile(outdir,'rcs_boresight_cuts.png'),'Resolution',220);
savefig(fig,fullfile(outdir,'rcs_boresight_cuts.fig'));
print(fig,fullfile(outdir,'rcs_boresight_cuts.svg'),'-dsvg'); close(fig);
fprintf('LOCAL CUTS: HH widths az %.6f deg, el %.6f deg.\n',azWidth.width_deg,elWidth.width_deg);
verify_pattern_peak;
end

function result=find_width(x,y,threshold)
x=x(:); y=y(:);
i=find(x<0 & y<=threshold,1,'last'); j=find(x>0 & y<=threshold,1,'first');
left=x(i)+(threshold-y(i))*(x(i+1)-x(i))/(y(i+1)-y(i));
right=x(j-1)+(threshold-y(j-1))*(x(j)-x(j-1))/(y(j)-y(j-1));
result=struct('left_deg',left,'right_deg',right,'width_deg',right-left,'threshold_dBsm',threshold);
end
