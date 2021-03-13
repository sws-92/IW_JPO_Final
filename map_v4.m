%% Plot Map
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
%%
clear
clear global

station=load('sinfit_temp_data.mat');
Tstn=load('/ocean/sstevens/IW_project/data/represent_stns_t.mat');
mod_lat=ncread('/ocean/eolson/MEOPAR/NEMO-forcing/grid/mesh_mask201702.nc',...
    'nav_lat'); mod_lat(mod_lat==0)=NaN;
mod_lon=ncread('/ocean/eolson/MEOPAR/NEMO-forcing/grid/mesh_mask201702.nc',...
    'nav_lon'); mod_lon(mod_lon==0)=NaN;
idx=contains(station.rel_station_names,{'GO-4';'PR-2'});
station.station_lon(idx)=[];station.station_lat(idx)=[];
station.rel_dataset(idx)=[];
thal=load('zage_plotting.mat','bottomY','bottom');
load NE_path.mat
DP_path=load('DP_path.mat');
fname='/ocean/rich/more/mmapbase/noaa_bc3/british_columbia_3_msl_2013.nc';

% load('BCcoast');

% This script plots the mean IW currents from the SSC model
load model_IW_outputUV.mat
load('BCcoast');
model_uv.lat(model_uv.lat==0)=NaN;
model_uv.lon(model_uv.lon==0)=NaN;

mean_u=mean(model_uv.u,3,'omitnan');
mean_v=mean(model_uv.v,3,'omitnan');

uv_mag=NaN(size(model_uv.u));

for i=1:size(model_uv.u,3)
        uv_mag(:,:,i)=(model_uv.u(:,:,i).^2+model_uv.v(:,:,i).^2).^0.5;
end
       
% Rotate data to create along/across strait velocities
% Conversion to from along/across to cartesian 
coastang=-30;
uv_i=mean_u+1i*mean_v;
mean_north=imag(uv_i.*exp(-1i*coastang*pi/180));
mean_east=real(uv_i.*exp(-1i*coastang*pi/180));

blockN=BlockMean(mean_north,15,15);blocklat=BlockMean(double(model_uv.lat),15,15);
blockE=BlockMean(mean_east,15,15);blocklon=BlockMean(double(model_uv.lon),15,15);

%%
figure('units','centimeters','outerposition',[0 0 19 21]);
lat_lim=[48.1 50.3]; lon_lim=[-125.99298287 -122.22345235];
ax1=axes;

% Set up figure
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');
caxis([-400 0]);
colm=m_colmap('blues');
colm=colm(15:end,:);
colormap(colm);
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');

ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
lnes=lines;
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

mp=m_patch([-123.6413;-123.5503;-122.3757;-122.8],...
    [48.5701;48.4196;48.4196;49.0969],rgb_x('light red'),'facealpha',0.25,...
    'linestyle','none');

m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none');
% m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','speckle','color',[1,1,1].*0.4);
% m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNWrivers.mat','linewidth',0.05,'color','w');
[Cmain,hmain]=m_contourf(lon(ilon),lat(ilat),Z'*-1,0:100:400,'linestyle','none');
caxis([0 400]);
cm=flipud(m_colmap('blues',7));cm([1 end-1:end],:)=[];
colormap(cm);

% Plot additional map information
m_line([-124.4331;-123.6],[49.1633;49.7210],'color','k','linestyle','--');
m_line([-123.6413;-122.8],[48.5701;49.0969],'color','k','linestyle','--');
m_line([-123.5503;-122.3757],[48.4196;48.4196],'color','k','linestyle','--');

m_text(-124.25,48.9,{'Vancouver';'Island'},'fontsize',8,...
    'horizontalalignment','center');
m_text(-124.4737,49.2828,'N.B.','fontsize',8,'Fontweight','bold');
m_text(-124.25,49.2164,'S.B.','fontsize',8,'Fontweight','bold');
m_text(-123.55,48.5074,'H.S.','fontsize',8,'Fontweight','bold');
m_text(-124.7483,48.5,'Juan de Fuca Strait','fontsize',8,'rotation',-24);
m_scatter(-123.1878,49.2518,35,'p','filled','markerfacecolor','k')
m_text(-123.1378,49.2518,'Vancouver','fontsize',8);
m_text(-124.475,49.7,'Texada Island','fontsize',6,'rotation',-40);
m_text(-124.3563,49.7506,'Malaspina Strait','fontsize',6,'rotation',-50);
m_text(-124.6793,49.7276,'Sabine Channel','fontsize',6,'rotation',-40);
m_text(-123.0446,49.0872,'Fraser River','fontsize',6,'rotation',30);
m_text(-125.3620,50.0768,{'Discovery';'Passage'},'fontsize',6,'rotation',-65);
m_text(-124.3336,49.5061,'L.I.','fontsize',6,'rotation',-30);
% m_annotation('textarrow',[-125.3783  -125.2400],[49.7903 49.9694],...
%     'string',{'Discovery';'Passage'},'fontsize',6,'Fontweight','normal',...
%     'horizontalalignment','center');
m_annotation('textarrow',[-122.6220  -123.0376],[48.9262 48.7513],...
    'string',{'Boundary';'Pass'},'fontsize',6,'Fontweight','normal',...
    'horizontalalignment','center');
% m_plot(NE_path.lon_list(1:20:end),NE_path.lat_list(1:20:end),...
%     'k--','linewidth',1);
% m_annotation('arrow',[NE_path.lon_list(end-20),NE_path.lon_list(end)],...
%     [NE_path.lat_list(end-20),NE_path.lat_list(end)],'color',rgb_x('black'),'linewidth',1);
% m_plot(mod_lon(DP_path.points(1:20:end)),mod_lat(DP_path.points(1:20:end)),...
%     'k--','linewidth',1);
% m_annotation('arrow',[mod_lon(DP_path.points(20)),mod_lon(DP_path.points(1))],...
%     [mod_lat(DP_path.points(20)),mod_lat(DP_path.points(1))],'color',rgb_x('black'),'linewidth',1);

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',8,'box','fancy');

drawnow;  
hFills = hmain.FacePrims;  
[hFills.ColorType] = deal('truecoloralpha'); 
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = 180;  
end

[CF,CFS]=m_contfbar(ax1,0.675,[0.8 0.95],Cmain,hmain,'endpiece','yes','axfrac',.03,'fontsize',8);
xlabel(CF,{'Depth (m)'},'fontsize',8,'color','k','fontweight','normal','rotation',0);
drawnow;  
hFills = CFS.FacePrims;  
[hFills.ColorType] = deal('truecoloralpha'); 
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = 180;  
end

ax2=axes('Position',[0.625 0.625 .4 0.35]);
lat_lim=[46.16487459 52];
lon_lim=[-127.5 -121];
m_proj('UTM','lon',lon_lim,'lat',lat_lim);  
[ELEV,ELON,ELAT]=m_etopo2([lon_lim lat_lim]);
ELEV(ELEV>0)=NaN;


fname2='/ocean/sstevens/IW_project/data/ubcSSn2DMeshMaskV17-02_f52b_5249_fe72.nc';
lat2=ncread(fname2,'gphit');
lon2=ncread(fname2,'glamt');
S=lat2==min(lat2(:));N=lat2==max(lat2(:));
W=lon2==min(lon2(:));E=lon2==max(lon2(:));
[lo,la]=m_ll2xy([lon2(E) lon2(N) lon2(W) lon2(S)]',...
    [lat2(E) lat2(N) lat2(W) lat2(S)]');
% [lo,la]=m_ll2xy(lon2,lat2);
as=alphaShape(lo(:),la(:));
[Cmain,hmain]=m_contourf(ELON,ELAT,ELEV*-1,0:100:400,'linestyle','none');
caxis([0 400]);
cm=flipud(m_colmap('blues',7));cm([1 end-1:end],:)=[];
colormap(cm);

hold on
% m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.6, 'Edgecolor', 'none');
m_gshhs_i('patch',rgb('light grey'),'edgecolor','none');
edges = boundaryFacets(as);
plot(lo(edges),la(edges),'k');
% plot(as,'FaceColor',rgb_x('light red'),'FaceAlpha',0.2,'linestyle','none');
m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
    'xaxisloc','top','yaxisloc','right','fontsize',7,'xtick',-127:2:-121);
m_rectangle(-125.99298287,48,3.7695,2.3,0,'edgecolor','r');
m_text(-122.3,49.2176,'Canada','fontsize',6');
m_text(-122.3,48.8346,'U.S.A.','fontsize',6');
m_text(-127,47,{'Pacific';'Ocean'},'fontsize',6');
m_text(-126.95,50.3,'J.S.','fontsize',6');
m_gshhs('fb1','linest','--','linewidth',1,'color',rgb_x('light red'));

drawnow;  
hFills = hmain.FacePrims;  
[hFills.ColorType] = deal('truecoloralpha'); 
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = 180;  
end
% 
% ax3=axes('Position',[0.01 0.025 0.4 0.4]);
% lat_lim=[48.4 50.15];
% lon_lim=[-125.280000012121 -122.500002324];
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
% hold on
% RGB=rgb('light grey');
% 
% [CS,CH]=m_contourf(model_uv.lon,model_uv.lat,...
%     mean_v,[-0.2:0.025:0.2],'linecolor','none');
% colormap(ax3,cmocean('balance'));
% % m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');
% hold on
% 
% idx=blockN>0;
% mv=m_vec(0.5,blocklon(idx),blocklat(idx),...
%     blockE(idx),blockN(idx),...
%     'shaftwidth',0.75,'headlength',5,'edgeclip', 'on','headangle',30,...
%     'facecolor',rgb_x('red'),'edgecolor','none');
% idx=blockN<0;
% mv=m_vec(0.5,blocklon(idx),blocklat(idx),...
%     blockE(idx),blockN(idx),...
%     'shaftwidth',0.75,'headlength',5,'edgeclip', 'on','headangle',30,...
%     'facecolor',rgb_x('blue'),'edgecolor','none');
% 
% m_gshhs_i('patch',rgb('light grey'),'edgecolor','none');
% 
% uv_arr=-0.2+1i*0;
% arr_along=imag(uv_arr.*exp(-1i*29*pi/180));
% arr_across=real(uv_arr.*exp(-1i*29*pi/180));
% m_vec(1,-122.9,49.95,arr_across,arr_along,...
%     'shaftwidth',0.75,'headlength',5,'edgeclip', 'on','headangle',30);
% m_text(-123.1,49.9,{'Along-strait direction';'(0.2 m s^{-1})'}','fontsize',7,...
%     'rotation',-29,'horizontalalignment','center');
% 
% m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','left','fontsize',6);
% [ax,h]=m_contfbar(-0.05,[0.125 0.875],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
% ylabel(ax,{'Along-Strait Current (m s ^{-1})'},'fontsize',8,'color','k','fontweight','bold');
% 
% ax4=axes('Position',[0.035 0.45 0.35 0.2]);
% cm=cmocean('grey',6);
% cm(1,:)=rgb_x('red');cm(2,:)=rgb_x('light red');
% cm(end-1,:)=rgb_x('light blue');
% cm(end,:)=[0.9 0.9 0.9];
% hold on
% ar(1)=rectangle('position',[thal.bottomY(end) 200 range(thal.bottomY) 250],...
%     'FaceColor',[0 0 0 0.55],'edgecolor','none');
% ar(2)=rectangle('position',[thal.bottomY(end) 50 49.5-thal.bottomY(end) 150],...
%     'FaceColor',cm(1,:),'edgecolor','none');
% ar(3)=rectangle('position',[49.5 50 thal.bottomY(1)-49.5 150],...
%     'FaceColor',cm(2,:),'edgecolor','none');
% ar(4)=rectangle('position',[thal.bottomY(end) 0 range(thal.bottomY) 50],...
%     'FaceColor',[cm(3,:) 0.55],'edgecolor','none');
% ar(5)=rectangle('position',[thal.bottomY(end) 0 0.1 450],...
%     'FaceColor',[cm(5,:) 1]);%,'edgecolor','none');
% ar2=area(thal.bottomY,thal.bottom,450,'FaceColor',rgb_x('black'));
% 
% % annotation('arrow',[0.81 0.76],[0.96 0.96]);
% % annotation('arrow',[0.81 0.76],[0.92 0.92]);
% set(gca,'ydir','reverse','YAxisLocation','right');
% axis tight;
% yticks([50 200])
% xticks([]);
%%
% 
% figure
% surf(ELON,ELAT,ELEV,'linestyle','none')
% 
% 
% idx=ELAT<48.5 | ELAT<49 & ELON<-124 | ELON<-125.5;
% tmpla=ELAT;tmpla(idx)=NaN;
% tmplo=ELON;tmplo(idx)=NaN;
% tmpE=ELEV;tmpE(idx)=NaN;
% 
% figure
% surf(tmplo,tmpla,tmpE,'linestyle','none')

%%
% set(gcf,'color','w');
% export_fig /ocean/sstevens/IW_project/figures/paper/map_v4.png -png -m3
% export_fig /ocean/sstevens/IW_project/figures/paper/vec/map_v4.eps -deps
set(gcf,'color','w');
export_fig /ocean/sstevens/IW_project/figures/paper/vec/map_v5.pdf -dpdf
% %%
% print /ocean/sstevens/IW_project/figures/paper/vec/map_v4.eps -depsc2 -painters
% 
% print /ocean/sstevens/IW_project/figures/paper/vec/map_v4_HQ.pdf -dpdf -painters
% print /ocean/sstevens/IW_project/figures/paper/vec/map_v4.pdf -dpdf
