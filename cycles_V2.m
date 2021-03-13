% Plot spatial variation in phase and amplitude
% addpath(genpath('/ocean/sstevens/'));
% addpath(genpath('/ocean/rich/home/matlab/m_map/'));

%%
clc
clear
clear global
surf=load('sinfit_temp_data_plotting.mat');
core=load('sinfit_temp_data_plotting_core.mat');
deep=load('sinfit_temp_data_plotting_deep.mat');
load('BCcoast');


lon_lim=[-127 -121];
lat_lim=[46.5 52];
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection

yl=47.85:0.01:51;
xl=ones(1,length(yl))*-124.0342;
[x,y]=m_ll2xy(xl,yl);
v=[x;y];
y_cent=y(round(length(y)/2));x_cent=x(round(length(y)/2));
R=[cosd(60) -sind(60);sind(60) cosd(60)];

AAS=[];
% Do rotation
[AAS(1,:),AAS(2,:)]=m_ll2xy(surf.station_lon,surf.station_lat);
center=repmat([x_cent; y_cent],1,length(surf.station_lon));
sa=AAS-center;so=R*sa;surf.AAS=(so+center)*6370;
surf.AAS(1,:)=surf.AAS(1,:)-min(surf.AAS(1,:));

AAS=[];
[AAS(1,:),AAS(2,:)]=m_ll2xy(core.station_lon,core.station_lat);
center=repmat([x_cent; y_cent],1,length(core.station_lon));
sa=AAS-center;so=R*sa;core.AAS=(so+center)*6370;
core.AAS(1,:)=core.AAS(1,:)-min(core.AAS(1,:));

AAS=[];
[AAS(1,:),AAS(2,:)]=m_ll2xy(deep.station_lon,deep.station_lat);
center=repmat([x_cent; y_cent],1,length(deep.station_lon));
sa=AAS-center;so=R*sa;deep.AAS=(so+center)*6370;
core.AAS(1,:)=core.AAS(1,:)-min(core.AAS(1,:));

% Lagrangian
load lagrangianM
% Find all ages in strait
load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;
inside_idx=inShape(station_bound,BP.lon_grid/100,BP.lat_grid/100);

straitA=BP.mean_map(inside_idx);
straitLa=BP.lat_grid(inside_idx)/100;
straitLo=BP.lat_grid(inside_idx)/100;

% create lines
center_line=[linspace(-125.15,-122.75,1000);linspace(50.15,48.4,1000)];
load centre_line_dist.mat

distance_idx=NaN(1,length(straitA));
Ldistance_N=NaN(1,length(straitA));
for i=1:length(straitA)
    [~,idxb]=min(abs(center_line(2,:)-straitLa(i)));
    Ldistance_N(i)=center_line_dist(idxb);
end

% dist=gsw_distance([ones(length(straitA),1)*-122.4364 straitLo],...
%     [ones(length(straitA),1)*48.1429 straitLa]);
% dist=(dist-min(dist))/1000;

Ldistance_N(isnan(straitA))=[];
straitA(isnan(straitA))=[];

[coeffs(1), coeffs(2)] = TheilSen([Ldistance_N'/1000,straitA]);
Lfittedy=polyval(coeffs,Ldistance_N/1000);

% model
m_core=load('sinfit_SSC_modelres_core.mat');
qlon=m_core.grid_lon;qlon(isnan(qlon))=0;
qlat=m_core.grid_lat;qlat(isnan(qlat))=0;

% Find all ages in strait
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;
inside_idx=inShape(station_bound,qlon,qlat);

straitA=m_core.cold_day(inside_idx);
straitLa=m_core.grid_lat(inside_idx);
msk=isnan(straitA); 
straitA(msk)=[];
straitLa(msk)=[];

% create lines
center_line=[linspace(-125.15,-122.75,1000);linspace(50.15,48.4,1000)];
load centre_line_dist.mat

distance_idx=NaN(1,length(straitA));
Mdistance_N=NaN(1,length(straitA));
for i=1:length(straitA)
    [~,idxb]=min(abs(center_line(2,:)-straitLa(i)));
    Mdistance_N(i)=center_line_dist(idxb);
end

[coeffs(1), coeffs(2)] = TheilSen([Mdistance_N'/1000,straitA]);
Mfittedy=polyval(coeffs,Mdistance_N/1000);

%% Plot seasonal cycles from all depths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot quick map
figure('units','centimeters','outerposition',[0 0 14 14],'color','w');
ax1=axes('Position',[0.1 0.11 .2 0.8]);
lat_lim=[50.15 48.4];
lon_lim=[-125.280000012121 -122.500002324];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.25) 
hold on
m_gshhs_f('patch',rgb('light grey'),'edgecolor','none');
% m_grid('linestyle','none','linewidth',1,'linecolor',rgb('light grey'),...
%     'tickdir','out','xaxisloc','bottom','yaxisloc','left','fontsize',8,...
%     'xtick',-126:-122,'ytick',48:51,'linealpha',0.05);
m_grid('linestyle',':','linewidth',1,'linecolor',rgb('very light grey'),...
    'tickdir','out','xaxisloc','top','yaxisloc','left','fontsize',8,...
    'xtick',-126:-122,'ytick',48:51,'linealpha',0.05);
text(0.85,0.1,'a)','units','normalized','fontweight','bold','fontsize',8);

% Plot stations and ruler
% m_scatter(surf.Dsort_lon,surf.Dsort_lat,25,'r','filled',...
%     'markerfacealpha',0.2,'markeredgealpha',0);
% m_plot(surf.center_line(1,:),surf.center_line(2,:),'k-');
    
% m_track(surf.center_line(1,:),surf.center_line(2,:),surf.center_line_dist/1000);

% cm_magma=magma(m);
% cm_inferno=inferno(m);
% cm_plasma=plasma(m);
% cm_viridis=viridis(m);

mycol=m_colmap('jet',length(surf.Dsort_coldday));
% mycol=flipud(viridis(length(surf.Dsort_coldday)+50));
% mycol(1:30,:)=[];

for i=1:length(surf.Dsort_coldday)
    m_scatter(surf.Dsort_lon(i),surf.Dsort_lat(i),25,mycol(i,:),'filled',...
        'markerfacealpha',0.6); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot seasonal cycles as contour plot

[X,Y]=meshgrid(surf.distance_N,1:365);
ax2=axes('Position',[0.4 0.645 .3 0.25]);
hold on
ax2.XAxis.Visible='off';
[CS,CH]=contourf(Y,X/1000,surf.Dsort_harms,[7:0.5:12],...
    'linecolor','none');
colormap(ax2,cmocean('balance'));
[ax,h]=m_contfbar(ax2,[0 1],1.1,CS,CH,'endpiece','no','axfrac',.1,'fontsize',8);
ylabel(ax,{'Temperature (^\circC)'},'fontsize',8,'color','k','fontweight','bold','rotation',0);
xticks(ax,7:12);
box off
% set(gca,'YColor','none')
mycol=m_colmap('jet',length(surf.Dsort_coldday));
for i=1:length(surf.Dsort_coldday)
    scatter(-10,surf.distance_N(i)/1000,30,mycol(i,:),'s','filled',...
        'markerfacealpha',0.5); 
end
xlim([-10 365]);
yticks(50:50:250);
text(0.85,0.1,'b)','units','normalized','fontweight','bold','fontsize',8);

ax3=axes('Position',[0.4 0.37 .3 0.25]);
[X,Y]=meshgrid(core.distance_N,1:365);
hold on
ax3.XAxis.Visible='off';
contourf(Y,X/1000,core.Dsort_harms,7:0.5:12,'linecolor','none');
colormap(ax3,cmocean('balance'));
box on
ylabel('Along-Strait Distance (km)','Fontsize',8,'fontweight','bold');
mycol=m_colmap('jet',length(core.Dsort_coldday));
for i=1:length(core.Dsort_coldday)
    scatter(-10,core.distance_N(i)/1000,30,mycol(i,:),'s','filled',...
        'markerfacealpha',0.5); 
end
yticks(50:50:250);
text(0.85,0.1,'e)','units','normalized','fontweight','bold','fontsize',8);

ax4=axes('Position',[0.4 0.095 .3 0.25]);
[X,Y]=meshgrid(deep.distance_N,1:365);
hold on
ax4.XAxis.Visible='off';
contourf(Y,X/1000,deep.Dsort_harms,7:0.5:12,'linecolor','none');
colormap(ax4,cmocean('balance'));
box on
mycol=m_colmap('jet',length(deep.Dsort_coldday));
for i=1:length(deep.Dsort_coldday)
    scatter(-10,deep.distance_N(i)/1000,30,mycol(i,:),'s','filled',...
        'markerfacealpha',0.5); 
end
yticks(50:50:250);
text(0.85,0.1,'h)','units','normalized','fontweight','bold','fontsize',8);

iax=axes('Position',[0.4 0.095 .3 0.8]);
linkaxes([iax ax2 ax3 ax4],'x'); 
iax.YAxis.Visible='off';
uistack(iax,'top');
iax.Color='none'; iax.XGrid='on';
% iax.XMinorGrid='on'; iax.MinorGridLineStyle='-';
iax.GridAlpha=0.05;
box on
xlabel('Yearday','fontsize',8,'fontweight','bold');
% iax.MinorGridAlpha=0.05;
xticks([0:100:300]);


%%% Add Coldday trend plots
% surface
axCD(1)=axes('Position',[0.72 0.645 .1 0.25]);
hold on

msk=surf.Dsort_cold_daySD<50;
errorbar(surf.Dsort_coldday(msk),surf.distance_N(msk)/1000,surf.Dsort_cold_daySD(msk)',...
    'horizontal','color',[0.8 0.8 0.8],'linewidth',1,'capsize',4);
scatter(surf.Dsort_coldday(msk),surf.distance_N(msk)/1000,5,'filled',...
    'markerfacecolor',[0.7 0.7 0.7]);

[coeffs(1),coeffs(2)]=tsreg(surf.distance_N/1000,surf.Dsort_coldday,...
    length(surf.Dsort_coldday));
fittedy=polyval(coeffs,surf.distance_N/1000);
err=sqrt((sum((fittedy-surf.distance_N).^2))/(length(fittedy)-2))*2/...
            sqrt(sum((surf.Dsort_coldday-nanmean(surf.Dsort_coldday)).^2));
plot(fittedy,surf.distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
axCD(1).YAxis.Visible='off';axCD(1).XAxis.Visible='off';
linkaxes([axCD(1) ax2],'y');
text(0.85,0.1,'c)','units','normalized','fontweight','bold','fontsize',8);

lnes=lines;
axCD(2)=axes('Position',[0.72 0.37 .1 0.25]);
hold on
msk=core.Dsort_cold_daySD<50;
errorbar(core.Dsort_coldday(msk),core.distance_N(msk)/1000,core.Dsort_cold_daySD(msk)',...
    'horizontal','color',[0.8 0.8 0.8],'linewidth',1,'capsize',4);
scatter(core.Dsort_coldday(msk),core.distance_N(msk)/1000,5,'filled',...
    'markerfacecolor',[0.7 0.7 0.7]);

[coeffs(1),coeffs(2)]=tsreg(core.distance_N/1000,core.Dsort_coldday,...
    length(core.Dsort_coldday));
fittedy=polyval(coeffs,core.distance_N/1000);
plot(Mfittedy,Mdistance_N/1000,'--','color',rgb('very light blue'),'linewidth',2);
plot(Lfittedy,Ldistance_N/1000,'--','color',rgb('very light orange'),'linewidth',2);
plot(fittedy,core.distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
axCD(2).YAxis.Visible='off';axCD(2).XAxis.Visible='off';
linkaxes([axCD(2) ax3],'y');
text(0.85,0.1,'f)','units','normalized','fontweight','bold','fontsize',8);

axCD(3)=axes('Position',[0.72 0.095 .1 0.25]);
hold on
msk=deep.Dsort_cold_daySD<50;
errorbar(deep.Dsort_coldday(msk),deep.distance_N(msk)/1000,deep.Dsort_cold_daySD(msk)',...
    'horizontal','color',[0.8 0.8 0.8],'linewidth',1,'capsize',4);
scatter(deep.Dsort_coldday(msk),deep.distance_N(msk)/1000,5,'filled',...
    'markerfacecolor',[0.7 0.7 0.7]);

[coeffs(1),coeffs(2)]=tsreg(deep.distance_N/1000,deep.Dsort_coldday,...
    length(deep.Dsort_coldday));
fittedy=polyval(coeffs,deep.distance_N/1000);
plot(fittedy,deep.distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
axCD(3).YAxis.Visible='off';axCD(3).XAxis.Visible='off';
linkaxes([axCD(3) ax4],'y');
text(0.85,0.1,'i)','units','normalized','fontweight','bold','fontsize',8);

iax2=axes('Position',[0.72 0.095 .1 0.8]);
linkaxes([iax2 axCD],'x'); 
iax2.YAxis.Visible='off';
uistack(iax2,'top');
iax2.Color='none'; iax2.XGrid='on';
iax2.GridAlpha=0.05;
xlabel('\phi_{CD} (days)','fontsize',8,'fontweight','bold');
% iax.MinorGridAlpha=0.05;
xticks([0:80:160]);
iax2.XMinorGrid='on'; iax2.MinorGridLineStyle='-';
iax2.XAxis.MinorTickValues = 0:40:200;
iax2.MinorGridAlpha=0.05;

%%% Add amp trend plots
% surface
axamp(1)=axes('Position',[0.85 0.645 .1 0.25]);
hold on

msk=surf.Dsort_ampSD<0.5;
errorbar(surf.Dsort_amp(msk),surf.distance_N(msk)/1000,surf.Dsort_ampSD(msk)',...
    'horizontal','color',[0.8 0.8 0.8],'linewidth',1,'capsize',4);
scatter(surf.Dsort_amp(msk),surf.distance_N(msk)/1000,5,'filled',...
    'markerfacecolor',[0.7 0.7 0.7]);

[coeffs(1),coeffs(2)]=tsreg(surf.distance_N/1000,surf.Dsort_amp,...
    length(surf.Dsort_amp));
fittedy=polyval(coeffs,surf.distance_N/1000);
plot(fittedy,surf.distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
axamp(1).YAxis.Visible='off';axamp(1).XAxis.Visible='off';
linkaxes([axamp(1) ax2],'y');
text(0.85,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);

axamp(2)=axes('Position',[0.85 0.37 .1 0.25]);
hold on

msk=core.Dsort_ampSD<0.5;
errorbar(core.Dsort_amp(msk),core.distance_N(msk)/1000,core.Dsort_ampSD(msk)',...
    'horizontal','color',[0.8 0.8 0.8],'linewidth',1,'capsize',4);
scatter(core.Dsort_amp(msk),core.distance_N(msk)/1000,5,'filled',...
    'markerfacecolor',[0.7 0.7 0.7]);

[coeffs(1),coeffs(2)]=tsreg(core.distance_N/1000,core.Dsort_amp,...
    length(core.Dsort_amp));
fittedy=polyval(coeffs,core.distance_N/1000);
plot(fittedy,core.distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
axamp(2).YAxis.Visible='off';axamp(2).XAxis.Visible='off';
linkaxes([axamp(2) ax3],'y');
text(0.85,0.1,'g)','units','normalized','fontweight','bold','fontsize',8);

axamp(3)=axes('Position',[0.85 0.095 .1 0.25]);
hold on

msk=deep.Dsort_ampSD<0.5;
errorbar(deep.Dsort_amp(msk),deep.distance_N(msk)/1000,deep.Dsort_ampSD(msk)',...
    'horizontal','color',[0.8 0.8 0.8],'linewidth',1,'capsize',4);
scatter(deep.Dsort_amp(msk),deep.distance_N(msk)/1000,5,'filled',...
    'markerfacecolor',[0.7 0.7 0.7]);

[coeffs(1),coeffs(2)]=tsreg(deep.distance_N/1000,deep.Dsort_amp,...
    length(deep.Dsort_amp));
fittedy=polyval(coeffs,deep.distance_N/1000);
plot(fittedy,deep.distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
axamp(3).YAxis.Visible='off';axamp(3).XAxis.Visible='off';
linkaxes([axamp(3) ax4],'y');
text(0.85,0.1,'j)','units','normalized','fontweight','bold','fontsize',8);

iax3=axes('Position',[0.85 0.095 .1 0.8]);
linkaxes([iax3 axamp],'x'); 
iax3.YAxis.Visible='off';
uistack(iax3,'top');
iax3.Color='none'; iax3.XGrid='on';
iax3.GridAlpha=0.05;
xlabel({'A_t (^\circC)'},'fontsize',8,'fontweight','bold');
% iax.MinorGridAlpha=0.05;
xticks([0:0.5:1.5]);
% iax2.XMinorGrid='on'; iax2.MinorGridLineStyle=':';
% iax2.XAxis.MinorTickValues = 0:40:200;
% iax2.MinorGridAlpha=0.05;

iax4=axes('Position',[0.4 0.645 .55 0.25]);
linkaxes([ax2 iax4],'y'); 
iax4.YAxis.Visible='off';iax4.XAxis.Visible='off';
uistack(iax4,'top');
iax4.Color='none'; iax4.YGrid='on';
iax4.GridAlpha=0.05;

iax5=axes('Position',[0.4 0.37 .55 0.25]);
linkaxes([ax3 iax5],'y'); 
iax5.YAxis.Visible='off';iax5.XAxis.Visible='off';
uistack(iax5,'top');
iax5.Color='none'; iax5.YGrid='on';
iax5.GridAlpha=0.05;

iax6=axes('Position',[0.4 0.095 .55 0.25]);
linkaxes([ax2 iax6],'y'); 
iax6.YAxis.Visible='off';iax6.XAxis.Visible='off';
uistack(iax6,'top');
iax6.Color='none'; iax6.YGrid='on';
iax6.GridAlpha=0.05;

set(findall(gcf,'-property','FontSize'),'FontSize',7)

%%
% export_fig /ocean/sstevens/IW_project/figures/paper/FITS_V5.png -png -m3
export_fig /ocean/sstevens/IW_project/figures/paper/vec/FITS_V7.pdf -dpdf
% print /ocean/sstevens/IW_project/figures/paper/vec/FITS_V5.pdf -dpdf -painters
%% Statistics
clc
[b1(1), b1(2)]=tsreg(surf.distance_N/1000,surf.Dsort_amp,...
    length(surf.Dsort_amp));
fittedy=polyval(b1,surf.distance_N/1000);
b1_err=sqrt((sum((fittedy-surf.Dsort_amp).^2))/(length(fittedy)-2))*2/...
    sqrt(sum((surf.distance_N/1000-nanmean(surf.distance_N/1000)).^2));

[b2(1), b2(2)]=tsreg(surf.Dsort_coldday,surf.distance_N,...
    length(surf.Dsort_coldday));
fittedy=polyval(b2,surf.distance_N);
b2_err=sqrt((sum((fittedy-surf.Dsort_coldday).^2))/(length(fittedy)-2))*2/...
    sqrt(sum((surf.distance_N-nanmean(surf.distance_N)).^2));

r1=corrcoef(surf.distance_N/1000,surf.Dsort_amp);
r2=corrcoef(surf.distance_N/1000,surf.Dsort_coldday);

fprintf('Surface amp=%2.3f+-%2.3f C 100km-1\nR=%2.2f\n\n',...
    b1(1)*100,b1_err*100,r1(2));
fprintf('Surface spd=%2.3f+-%2.3f ms-1\n R=%2.2f\n\n',...
    b2(1)/86000,b2_err/86000,r2(2)); 


[b1(1), b1(2)]=tsreg(core.distance_N/1000,core.Dsort_amp,...
    length(core.Dsort_amp));
fittedy=polyval(b1,core.distance_N/1000);
b1_err=sqrt((sum((fittedy-core.Dsort_amp).^2))/(length(fittedy)-2))*2/...
    sqrt(sum((core.distance_N/1000-nanmean(core.distance_N/1000)).^2));

[b2(1), b2(2)]=tsreg(core.Dsort_coldday,core.distance_N,...
    length(core.Dsort_coldday));
fittedy=polyval(b2,core.distance_N);
b2_err=sqrt((sum((fittedy-core.Dsort_coldday).^2))/(length(fittedy)-2))*2/...
    sqrt(sum((core.distance_N-nanmean(core.distance_N)).^2));

r1=corrcoef(core.distance_N/1000,core.Dsort_amp);
r2=corrcoef(core.distance_N/1000,core.Dsort_coldday);

fprintf('Core amp=%2.3f+-%2.3f C 100km-1\n R=%2.2f\n\n',...
    b1(1)*100,b1_err*100,r1(2));
fprintf('Core spd=%2.3f+-%2.3f ms-1\n R=%2.2f\n\n',...
    b2(1)/86000,b2_err/86000,r2(2)); 


[b1(1), b1(2)]=tsreg(deep.distance_N/1000,deep.Dsort_amp,...
    length(deep.Dsort_amp));
fittedy=polyval(b1,deep.distance_N/1000);
b1_err=sqrt((sum((fittedy-deep.Dsort_amp).^2))/(length(fittedy)-2))*2/...
    sqrt(sum((deep.distance_N/1000-nanmean(deep.distance_N/1000)).^2));

[b2(1), b2(2)]=tsreg(deep.Dsort_coldday,deep.distance_N,...
    length(deep.Dsort_coldday));
fittedy=polyval(b2,deep.distance_N);
b2_err=sqrt((sum((fittedy-deep.Dsort_coldday).^2))/(length(fittedy)-2))*2/...
    sqrt(sum((deep.distance_N-nanmean(deep.distance_N)).^2));

r1=corrcoef(deep.distance_N/1000,deep.Dsort_amp);
r2=corrcoef(deep.distance_N/1000,deep.Dsort_coldday);

fprintf('Deep amp=%2.3f+-%2.3f C 100km-1\n R=%2.2f\n\n',...
    b1(1)*100,b1_err*100,r1(2));
fprintf('Deep spd=%2.3f+-%2.3f ms-1\n R=%2.2f\n\n',...
    b2(1)/86000,b2_err/86000,r2(2)); 

% Amplitude attenuation
amp_low=min(surf.Dsort_amp(surf.distance_N/1000>200));
amp_high=min(surf.Dsort_amp(surf.distance_N/1000<50));

p(1)=100-(amp_low/amp_high)*100;

amp_low=min(core.Dsort_amp(core.distance_N/1000>200));
amp_high=min(core.Dsort_amp(core.distance_N/1000<50));

p(2)=100-(amp_low/amp_high)*100;

amp_low=min(deep.Dsort_amp(deep.distance_N/1000>200));
amp_high=min(deep.Dsort_amp(deep.distance_N/1000<50));

p(3)=100-(amp_low/amp_high)*100;

fprintf('Amplitude atten. pct.=%2.0f-%2.0f\n\n',...
    min(p),max(p)); 
