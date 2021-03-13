%% Supplementary Information- IW paper
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
%%
clear
load /ocean/sstevens/IW_project/data/sinfit_ox.mat
%
figure('units','centimeters','outerposition',[0 0 20 10],'color','w');
ax1=axes('Position',[0.0 0.11 .2 0.8]);
lat_lim=[50.15 48.4];
lon_lim=[-125.280000012121 -122.500002324];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.25) 
hold on
m_gshhs_f('patch',rgb('light grey'),'edgecolor','none');
m_grid('linestyle',':','linewidth',1,'linecolor',rgb('very light grey'),...
    'tickdir','out','xaxisloc','top','yaxisloc','left','fontsize',8,...
    'xtick',-126:-122,'ytick',48:51,'linealpha',0.05);
text(0.75,0.1,'a)','units','normalized','fontweight','bold','fontsize',8);
mycol=m_colmap('jet',length(Dsort_min_day));

for i=1:length(Dsort_min_day)
    m_scatter(Dsort_lon(i),Dsort_lat(i),25,mycol(i,:),'filled',...
        'markerfacealpha',0.6); 
end

axes('Position',[0.25 0.11 .15 0.8]);
[X,Y]=meshgrid(distance_N,1:365);
hold on
[CS,CH]=contourf(Y,X/1000,Dsort_harms,...
    'linecolor','none');
colormap(cmocean('balance'));
[ax,h]=m_contfbar(1.1,[0 1],CS,CH,'endpiece','no','axfrac',.1,'fontsize',8);
ylabel(ax,{'D.O. (\mumol kg^{-1})'},'fontsize',8,'color','k','fontweight','bold','rotation',-90);
xticks(ax,7:12);
box off
mycol=m_colmap('jet',length(Dsort_min_day));
for i=1:length(Dsort_min_day)
    scatter(-10,distance_N(i)/1000,60,mycol(i,:),'s','filled',...
        'markerfacealpha',0.5); 
end
xlim([-10 365]);
yticks(50:50:250);
xlabel('Yearday','fontsize',8,'fontweight','bold');
ylabel('Along-Strait Excursion (km)','Fontsize',8,'fontweight','bold');
text(0.85,0.1,'b)','units','normalized','fontweight','bold','fontsize',8);

axes('Position',[0.55 0.11 .15 0.35]);
scatter(Dsort_min_day,distance_N/1000,5,'filled',...
    'markerfacecolor',[0.8 0.8 0.8]);
hold on
[coeffs(1),coeffs(2)]=tsreg(distance_N/1000,Dsort_min_day,...
    length(Dsort_min_day));
fittedy=polyval(coeffs,distance_N/1000);
plot(fittedy,distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
ylabel('Along-Strait Excursion (km)','fontsize',7,'fontweight','bold');
xlabel('Least oxic yearday','fontsize',7,'fontweight','bold');
text(0.85,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);

axes('Position',[0.55 0.55 .15 0.35]);
scatter(Dsort_oamp,distance_N/1000,5,'filled',...
    'markerfacecolor',[0.8 0.8 0.8]);
hold on
[coeffs(1),coeffs(2)]=tsreg(distance_N/1000,Dsort_oamp,...
    length(Dsort_min_day));
fittedy=polyval(coeffs,distance_N/1000);
plot(fittedy,distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
xlabel('Amplitude (\mumol kg^{-1})','fontsize',7,'fontweight','bold');
text(0.85,0.1,'c)','units','normalized','fontweight','bold','fontsize',8);

axes('Position',[0.8 0.11 .15 0.8]);
scatter(Dsort_oxmean,distance_N/1000,5,'filled',...
    'markerfacecolor',[0.8 0.8 0.8]);
hold on
[coeffs(1),coeffs(2)]=tsreg(distance_N/1000,Dsort_oxmean,...
    length(Dsort_min_day));
fittedy=polyval(coeffs,distance_N/1000);
plot(fittedy,distance_N/1000,'color',rgb('light red'),'linewidth',2);
axis tight
ylabel('Along-Strait Excursion (km)','fontsize',7,'fontweight','bold');
xlabel('Mean D.O. (\mumol kg^{-1})','fontsize',7,'fontweight','bold');
text(1,0.1,'e)','units','normalized','fontweight','bold','fontsize',8);

set(findall(gcf,'-property','FontSize'),'FontSize',8)

%%
export_fig /ocean/sstevens/IW_project/figures/paper/SI/ox_fits.png -m3

% Find DO utilization
[coeffs(1),coeffs(2)]=tsreg(Dsort_min_day,Dsort_oxmean,...
    length(Dsort_min_day));
fittedy=polyval(coeffs,Dsort_min_day);
err=sqrt((sum((fittedy-Dsort_oxmean).^2))/(length(fittedy)-2))*2/...
    sqrt(sum((Dsort_min_day-nanmean(Dsort_min_day)).^2));

%% set up kriging of coldest yearday
load /ocean/sstevens/IW_project/data/thalweg.mat
load('BCcoast');

% create grids and lines for kriging
west_line=linspace(-125.9,-123.5,1000);
east_line=linspace(-124.4,-122,1000);
center_line=[linspace(-125.15,-122.75,1000);linspace(50.3,48.5,1000)];

% center_line_dist=NaN(1,1000);
% for i=1:length(center_line)
%     center_line_dist(i)=gsw_distance([center_line(1,i) center_line(1,end)],...
%         [center_line(2,i) center_line(2,end)]);
% end
load centre_line_dist.mat

grid_lat=repmat([linspace(50.3,48.5,1000)]',1,1000);
grid_lon=NaN(1,1000);

% Create slanted grid
for i=1:length(grid_lat)
    grid_lon(i,:)=linspace(west_line(i),east_line(i),1000);
end

clc
fprintf('The grid has %3.2f m by %3.2f m spacing (lon x lat)',...
    gsw_distance([grid_lon(1,1),grid_lon(1,2)],...
    [grid_lat(1,1),grid_lat(1,1)]),...
    gsw_distance([grid_lon(1,1),grid_lon(1,1)],...
    [grid_lat(1,1),grid_lat(2,1)]));

lat_lim=[48.5 50.3];
lon_lim=[-125.280000012121 -123.000002324];

% Create shallow water mask for use in Kriging
fname='/ocean/rich/more/mmapbase/noaa_bc3/british_columbia_3_msl_2013.nc';
Zlat=ncread(fname,'lat');
Zlon=ncread(fname,'lon');
ilon=Zlon>=lon_lim(1) & Zlon<=lon_lim(2);
ilat=Zlat>=lat_lim(1) & Zlat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);
Z=fliplr(rot90(Z,3));
[Zlon_grid,Zlat_grid]=meshgrid(Zlon(ilon),Zlat(ilat));
Zint=interp2(Zlon_grid,Zlat_grid,Z,grid_lon,grid_lat);

% Set up XY projection
lat_limV=[min(grid_lat(:)) max(grid_lat(:))];
lon_limV=[min(grid_lon(:)) max(grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(station_lon,station_lat);
[Xgrid,Ygrid]=m_ll2xy(grid_lon,grid_lat);
RGB=rgb('light grey');

% Plotting
% dday = variogram([X' Y'],min_day','nrbins',50,'maxdist',1e5,'plotit',false);
dday = variogram([X' Y'],min_day');
[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);
zi = kriging_use(vstruct,station_lon',station_lat',min_day',grid_lon,grid_lat);

% Find youngest water in Haro Strait
HS_AS=load('HS_AS.mat');
AS=alphaShape(HS_AS.l,HS_AS.ll);
msk=inShape(AS,station_lon,station_lat);
zi=zi-min(min_day(msk));

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25],[station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

lat_lim=[50.42283746 48.252093478];
lon_lim=[-125.82436392 -122.252345235];

% figure('units','centimeters','outerposition',[0 0 19 23],'color','w'); % paper
figure('units','centimeters','outerposition',[0 0 17 22],'color','w');% poster
ax1=axes('Position',[0.035 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
caxis([0 120]);
[CS,CH]=m_contourf(grid_lon,grid_lat,zi,0:10:120,...
        'linestyle','none','linecolor','k','linewidth',0.1);
colormap(ax1,turbo);

m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none'); 
s=m_scatter(station_lon,station_lat,1,'x','k');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',7);

text(0.9,0.025,'a)','units','normalized','fontsize',8,'Fontweight','bold');
% text(0.05,0.95,'Age','units','normalized','fontsize',8,'Fontweight','bold');
% m_plot(pth.Xlon,pth.Xlat,'w:','linewidth',2);
% m_text(-123.0690,49.0852,sprintf('%2.1f cm s^{-1}',surf.EBCspeed),...
%     'fontsize',6,'color','w','fontweight','bold');
m_contour(grid_lon,grid_lat,Zint,[-10 -10],RGB,'linewidth',0.1);
[ax,h]=m_contfbar(ax1,[0.1 0.9],1.05,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',0:20:120,'linestyle','none');
xlabel(ax,'Age (days)','fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dday = variogram([X' Y'],oamp');
[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);
zi = kriging_use(vstruct,station_lon',station_lat',oamp',grid_lon,grid_lat);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25],[station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

lat_lim=[50.42283746 48.252093478];
lon_lim=[-125.82436392 -122.252345235];

% figure('units','centimeters','outerposition',[0 0 19 23],'color','w'); % paper
ax2=axes('Position',[0.4 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
caxis([10 60]);
[CS,CH]=m_contourf(grid_lon,grid_lat,zi,20:2.5:60,...
        'linestyle','none','linecolor','k','linewidth',0.1);
colormap(ax2,flipud(turbo));

m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none'); 
s=m_scatter(station_lon,station_lat,1,'x','k');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',7);

text(0.9,0.025,'b)','units','normalized','fontsize',8,'Fontweight','bold');
% text(0.05,0.95,'Age','units','normalized','fontsize',8,'Fontweight','bold');
% m_plot(pth.Xlon,pth.Xlat,'w:','linewidth',2);
% m_text(-123.0690,49.0852,sprintf('%2.1f cm s^{-1}',surf.EBCspeed),...
%     'fontsize',6,'color','w','fontweight','bold');
m_contour(grid_lon,grid_lat,Zint,[-10 -10],RGB,'linewidth',0.1);
[ax,h]=m_contfbar(ax2,[0.1 0.9],1.05,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',10:10:60,'linestyle','none');
xlabel(ax,'Amplitude (\mumol kg^{-1})','fontweight','bold');
%%
export_fig /ocean/sstevens/IW_project/figures/paper/SI/ox_maps.png -m3

%% SSC 
clear
load sinfit_SSC_plottting_LAYER.mat
%
clc
fprintf('The grid has %3.2f m by %3.2f m spacing (lon x lat)',...
    gsw_distance([grid_lon(1,1),grid_lon(1,2)],...
    [grid_lat(1,1),grid_lat(1,1)]),...
    gsw_distance([grid_lon(1,1),grid_lon(1,1)],...
    [grid_lat(1,1),grid_lat(2,1)]));

% Set up XY projection
lat_limV=[min(grid_lat(:)) max(grid_lat(:))];
lon_limV=[min(grid_lon(:)) max(grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(obs.station_lon,obs.station_lat);
[Xgrid,Ygrid]=m_ll2xy(grid_lon,grid_lat);

% dday = variogram([X' Y'],min_day','nrbins',50,'maxdist',1e5,'plotit',false);
dday = variogram([X' Y'],cold_day');
[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);
zi = kriging_use(vstruct,obs.station_lon',obs.station_lat',cold_day',grid_lon,grid_lat);

% Find youngest water in Haro Strait
HS_AS=load('HS_AS.mat');
AS=alphaShape(HS_AS.l,HS_AS.ll);
msk=inShape(AS,obs.station_lon,obs.station_lat);
zi=zi-min(cold_day(msk));

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([obs.station_lon';-125.25],[obs.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK
RGB=rgb('light grey');

lat_lim=[50.42283746 48.252093478];
lon_lim=[-125.82436392 -122.252345235];

% figure('units','centimeters','outerposition',[0 0 19 23],'color','w'); % paper
figure('units','centimeters','outerposition',[0 0 17 22],'color','w');% poster
ax1=axes('Position',[0.035 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
caxis([0 140]);
[CS,CH]=m_contourf(grid_lon,grid_lat,zi,0:10:140,...
        'linestyle','none','linecolor','k','linewidth',0.1);
colormap(ax1,mod_blue_orange);

m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none'); 
s=m_scatter(obs.station_lon,obs.station_lat,1,'x','k');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',7);

text(0.9,0.025,'a)','units','normalized','fontsize',8,'Fontweight','bold');
% text(0.05,0.95,'Age','units','normalized','fontsize',8,'Fontweight','bold');
% m_plot(pth.Xlon,pth.Xlat,'w:','linewidth',2);
% m_text(-123.0690,49.0852,sprintf('%2.1f cm s^{-1}',surf.EBCspeed),...
%     'fontsize',6,'color','w','fontweight','bold');
m_contour(grid_lon,grid_lat,Zint,[-10 -10],RGB,'linewidth',0.1);
[ax,h]=m_contfbar(ax1,[0.1 0.9],1.05,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',0:20:140,'linestyle','none');
xlabel(ax,'Age (days)','fontweight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cm1=cmocean('balance');
cm2=cmocean('curl');
cm=[cm1(20:110,:);cm2(140:end,:)];

dday = variogram([X' Y'],amp');
[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);
zi = kriging_use(vstruct,obs.station_lon',obs.station_lat',amp',grid_lon,grid_lat);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([obs.station_lon';-125.25],[obs.station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

lat_lim=[50.42283746 48.252093478];
lon_lim=[-125.82436392 -122.252345235];

% figure('units','centimeters','outerposition',[0 0 19 23],'color','w'); % paper
ax2=axes('Position',[0.4 0.05 0.3 0.8]);
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
hold on
caxis([0.4 1.6]);
[CS,CH]=m_contourf(grid_lon,grid_lat,zi,0:0.05:1.4,...
        'linestyle','none','linecolor','k','linewidth',0.1);
colormap(ax2,cm);

m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none'); 
s=m_scatter(obs.station_lon,obs.station_lat,1,'x','k');
m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',7);

text(0.9,0.025,'b)','units','normalized','fontsize',8,'Fontweight','bold');
% text(0.05,0.95,'Age','units','normalized','fontsize',8,'Fontweight','bold');
% m_plot(pth.Xlon,pth.Xlat,'w:','linewidth',2);
% m_text(-123.0690,49.0852,sprintf('%2.1f cm s^{-1}',surf.EBCspeed),...
%     'fontsize',6,'color','w','fontweight','bold');
m_contour(grid_lon,grid_lat,Zint,[-10 -10],RGB,'linewidth',0.1);
[ax,h]=m_contfbar(ax2,[0.1 0.9],1.05,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',0:0.2:1.4,'linestyle','none');
xlabel(ax,'Amplitude (^\circC)','fontweight','bold');
%%
export_fig /ocean/sstevens/IW_project/figures/paper/SI/SSC_station_maps.png -m3
