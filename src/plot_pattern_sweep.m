function plot_pattern_sweep(outdir)
% Redraw from saved actual SBR results; never substitutes an analytical map.
if nargin==0, outdir=pattern_output_dir; end
D=load(fullfile(outdir,'rcs_pattern_data.mat')); S=D.summary;
set(groot,'defaultFigureVisible','off','defaultAxesFontName','Microsoft YaHei',...
    'defaultTextFontName','Microsoft YaHei','defaultAxesFontSize',12);
[AZ,EL]=meshgrid(D.az,D.el);
HH=10*log10(D.HH); VV=10*log10(D.VV);
floorDB=-25; topDB=13; colorLimits=[floorDB,topDB];
fig=figure('Color','w','Position',[50,50,1250,570]);
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for pol=1:2
    ax=nexttile; if pol==1, Z=HH; name='HH'; else, Z=VV; name='VV'; end
    imagesc(D.az,D.el,Z); axis xy; axis square; hold on;
    contour(AZ,EL,Z,[S.boresight_HH_dBsm-3,S.boresight_HH_dBsm-3],...
        'k-','LineWidth',1.25);
    plot(45,S.boresight_el_deg,'wo','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',6);
    text(47,S.boresight_el_deg+3,sprintf('%.2f dBsm',S.boresight_HH_dBsm),'FontSize',11,'Color','k');
    xlim([0,90]); ylim([0,90]); xticks(0:15:90); yticks(0:15:90);
    xlabel('方位角 (°)'); ylabel('俯仰角 (°)'); title([name ' 同极化']);
    caxis(colorLimits); colormap(ax,parula(256)); box on;
end
cb=colorbar; cb.Layout.Tile='east'; cb.Label.String='单站 RCS (dBsm)';
title(tl,{'340 GHz 三面角 RCS：内棱长 43.6 mm，PEC SBR',...
    '第一卦限 1°–89°，步长 0.5°；黑线为相对正向值 −3 dB；色标下限 −25 dBsm'},...
    'FontSize',14,'FontWeight','normal');
exportgraphics(fig,fullfile(outdir,'rcs_heatmaps.png'),'Resolution',220);
savefig(fig,fullfile(outdir,'rcs_heatmaps.fig')); close(fig);

fig=figure('Color','w','Position',[50,50,1250,570]);
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
colors=[0.05,0.35,0.7;0.82,0.32,0.08];
for plane=1:2
    nexttile;
    if plane==1
        x=D.cutAz-45; y=10*log10(D.cutAzValues); xlab='方位角偏移 (°)，俯仰角固定 35.264°';
        w=S.HH_azimuth_3dB;
        heading=sprintf('方位切面：HH −3 dB 宽度 %.2f°',w.width_deg);
    else
        x=D.cutEl-S.boresight_el_deg; y=10*log10(D.cutElValues);
        xlab='俯仰角偏移 (°)，方位角固定 45°'; w=S.HH_elevation_3dB;
        heading=sprintf('俯仰切面：HH −3 dB 宽度 %.2f°',w.width_deg);
    end
    plot(x,y(:,1),'-','LineWidth',1.6,'Color',colors(1,:)); hold on;
    plot(x,y(:,2),'--','LineWidth',1.4,'Color',colors(2,:));
    yline(S.boresight_HH_dBsm-3,':','−3 dB','LabelHorizontalAlignment','left','LineWidth',1.2);
    xline(0,':','LineWidth',1); grid on; box on;
    xlabel(xlab); ylabel('单站 RCS (dBsm)'); title(heading);
    xlim([min(x),max(x)]); ylim([-35,15]);
    legend('HH','VV','Location','south');
end
title(tl,{'340 GHz 正向切面：正向为方位角 45°、俯仰角 35.264°',...
    '角度步长 0.25°；每个角度 589,824 个孔口采样点；未加入显式 UTD/PTD 边缘绕射'},...
    'FontSize',14,'FontWeight','normal');
exportgraphics(fig,fullfile(outdir,'rcs_principal_cuts.png'),'Resolution',220);
savefig(fig,fullfile(outdir,'rcs_principal_cuts.fig'));
print(fig,fullfile(outdir,'rcs_principal_cuts.svg'),'-dsvg'); close(fig);

fig=figure('Color','w','Position',[50,50,850,760]);
radius=sqrt(D.HH/max(D.HH(:)));
X=radius.*cosd(EL).*cosd(AZ); Y=radius.*cosd(EL).*sind(AZ); Z=radius.*sind(EL);
surf(X,Y,Z,HH,'EdgeColor','none','FaceColor','interp'); hold on;
n0=[1,1,1]/sqrt(3); plot3([0,n0(1)*1.04],[0,n0(2)*1.04],[0,n0(3)*1.04],'k--','LineWidth',1.3);
plot3(n0(1),n0(2),n0(3),'ko','MarkerFaceColor','w','MarkerSize',7);
axis equal; xlim([0,1.05]); ylim([0,1.05]); zlim([0,1.05]);
xticks(0:0.25:1); yticks(0:0.25:1); zticks(0:0.25:1);
view(125,25); grid on; box on; colormap(parula(256)); caxis(colorLimits);
cb=colorbar; cb.Label.String='HH 单站 RCS (dBsm)';
xlabel('x（归一化）'); ylabel('y（归一化）'); zlabel('z（归一化）');
title({'340 GHz 三面角 HH 三维方向图（仅第一卦限）',...
    '半径 = √(σ / σ_{max})；颜色 = 绝对 RCS；虚线为正向'},...
    'FontSize',14,'FontWeight','normal');
exportgraphics(fig,fullfile(outdir,'rcs_3d_HH.png'),'Resolution',220);
savefig(fig,fullfile(outdir,'rcs_3d_HH.fig')); close(fig);
end
