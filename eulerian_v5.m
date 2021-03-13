%% Plot Eularian Map- SSC, nodes, timescales
%addpath(genpath('/ocean/sstevens/'));
%addpath(genpath('/ocean/rich/home/matlab/m_map/'));
%%
clear
mod_lat=ncread('/ocean/sallen/allen/research/MEOPAR/grid/mesh_mask201702.nc',...
    'nav_lat'); mod_lat(mod_lat==0)=NaN;
mod_lon=ncread('/ocean/sallen/allen/research/MEOPAR/grid/mesh_mask201702.nc',...
    'nav_lon'); mod_lon(mod_lon==0)=NaN;
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

blockN=BlockMean(mean_north,10,10);blocklat=BlockMean(double(model_uv.lat),10,10);
blockE=BlockMean(mean_east,10,10);blocklon=BlockMean(double(model_uv.lon),10,10);

stacey=load('/ocean/rich/home/metro/currents/StaceyMean.mat');

%% Overall current comparison
mlat=ncread('/ocean/sallen/allen/research/MEOPAR/grid/mesh_mask201702.nc',...
    'nav_lat');
mlon=ncread('/ocean/sallen/allen/research/MEOPAR/grid/mesh_mask201702.nc',...
    'nav_lon');

ss=load('SoG_AS.mat');
AS=alphaShape(ss.x,ss.y);
AS.Alpha=AS.Alpha*2;

up=nanmean(mean_v(inShape(AS,mlon,mlat) & mean_v>0),'all');
down=nanmean(mean_v(inShape(AS,mlon,mlat) & mean_v<0),'all');

upM=max(mean_v(inShape(AS,mlon,mlat) & mean_v>0));
downM=min(mean_v(inShape(AS,mlon,mlat) & mean_v<0));

prop=up-abs(down)

%% Load in nodes
central=load('/ocean/kstankov/ADCP/central/ADCP_central.mat');
idx=central.chartdepth>50 & central.chartdepth<150;
central.mnV=nanmean(central.vtrue(idx,:),'all')/100;
central.mnU=nanmean(central.utrue(idx,:),'all')/100;
central.mnM=sqrt(central.mnV^2+central.mnU^2);
central.lat=49.0404366;
central.lon=-123.425795;
Idx=dsearchn([mod_lon(:) mod_lat(:)],[central.lon central.lat]);
col=zeros(length(mod_lat(:)),1);col(Idx)=1;
igrid=logical(reshape(col,size(mod_lat)));
[I,J]=find(igrid==1);
central.modU=mean_east(I,J);
central.modV=mean_north(I,J);

DDL=load('/ocean/kstankov/ADCP/ddl/ADCP_ddl.mat');
idx=DDL.chartdepth>50 & DDL.chartdepth<150;
DDL.mnV=nanmean(DDL.vtrue(idx,:),'all')/100;
DDL.mnU=nanmean(DDL.utrue(idx,:),'all')/100;
DDL.mnM=sqrt(DDL.mnV^2+DDL.mnU^2);
DDL.lat=49.085095;
DDL.lon=-123.330155;
Idx=dsearchn([mod_lon(:) mod_lat(:)],[DDL.lon DDL.lat]);
col=zeros(length(mod_lat(:)),1);col(Idx)=1;
igrid=logical(reshape(col,size(mod_lat)));
[I,J]=find(igrid==1);
DDL.modU=mean_east(I,J);
DDL.modV=mean_north(I,J);

east=load('/ocean/kstankov/ADCP/east/ADCP_east.mat');
idx=east.chartdepth>50 & east.chartdepth<150;
east.mnV=nanmean(east.vtrue(idx,:),'all')/100;
east.mnU=nanmean(east.utrue(idx,:),'all')/100;
east.mnM=sqrt(east.mnV^2+east.mnU^2);
east.lat=49.042835;
east.lon=-123.317265;
Idx=dsearchn([mod_lon(:) mod_lat(:)],[east.lon east.lat]);
col=zeros(length(mod_lat(:)),1);col(Idx)=1;
igrid=logical(reshape(col,size(mod_lat)));
[I,J]=find(igrid==1);
east.modU=mean_east(I,J);
east.modV=mean_north(I,J);

for i=1:length(stacey.C100)
    Idx=dsearchn([mod_lon(:) mod_lat(:)],[stacey.mLon(i) stacey.mLat(i)]);
    col=zeros(length(mod_lat(:)),1);col(Idx)=1;
    igrid=logical(reshape(col,size(mod_lat)));
    [I,J]=find(igrid==1);
    stacey.modC100U(i)=mean_east(I,J);
    stacey.modC100V(i)=mean_north(I,J);
end

stacey.mod_mag=sqrt(stacey.modC100U.^2+stacey.modC100V.^2);

save '/ocean/sstevens/IW_project/data/eulerianM.mat' blockE blockN blockE ...
     blockN mean_east mean_north mean_u mean_v mod_lon mod_lat
%%  Make Eulerian mean map
figure('units','centimeters','outerposition',[0 0 14 18],'color','w');% lat_lim=[48.4 50.15];
ax1=axes('Position',[0.05 0.55 0.9 0.45]);
mean_v(mean_v==0)=NaN;

cm=m_colmap('diverging',256);ct=round(length(cm)/2);
cm(ct-30:ct+30,:)=[];

lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -122.500002324];
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

caxis([-0.2 0.2]);
[CS,CH]=m_contourf(model_uv.lon,model_uv.lat,...
    mean_v,'linecolor','none','levelstep',0.01);
% colormap(ax1,cmocean('balance','pivot',0));
colormap(ax1,cm);
[ax,h]=m_contfbar(ax1,0.15,[0.525 0.9],CS,[-0.2:0.01:0.2],'endpiece','yes','axfrac',.04,...
    'fontsize',7,'edgecolor','none','levels','match');
% m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');
hold on

idx=blockN>0;
m_vec(0.5,blocklon(idx),blocklat(idx),...
    blockE(idx),blockN(idx),...
    'shaftwidth',0.75,'headlength',5,'edgeclip', 'on','headangle',30,...
    'facecolor',rgb_x('red'),'edgecolor','none');
idx=blockN<0;
m_vec(0.5,blocklon(idx),blocklat(idx),...
    blockE(idx),blockN(idx),...
    'shaftwidth',0.75,'headlength',5,'edgeclip', 'on','headangle',30,...
    'facecolor',rgb_x('blue'),'edgecolor','none');

m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none');
%m_gshhs_l('patch',rgb('light grey'),'edgecolor','none');

uv_arr=-0.2+1i*0;
arr_along=imag(uv_arr.*exp(-1i*29*pi/180));
arr_across=real(uv_arr.*exp(-1i*29*pi/180));
m_vec(0.5,-122.825,49.825,arr_across,arr_along,...
    'shaftwidth',0.75,'headlength',5,'edgeclip', 'on','headangle',30);
m_text(-123.05,49.75,{'Along-strait direction';'(0.2 m s^{-1})'}','fontsize',8,...
    'rotation',-29,'horizontalalignment','center');

m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',7,'xticklabels',[]);

ylabel(ax,{'Along-Strait';'Current (m s^{-1})'},'fontsize',7,'color','k','fontweight','bold');
m_rectangle(-123.6344,49,.6,.3,1,'color','k');
text(0.925,0.95,'a)','units','normalized','fontsize',8,'Fontweight','bold');

ax2=axes('Position',[0.025 0.45 .4 0.35]);
lat_lim=[49 49+0.275];
lon_lim=[-123.7 -123.7+0.6];
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on

[CS,CH]=m_contourf(model_uv.lon,model_uv.lat,...
    mean_v,[-0.2:0.01:0.2],'linecolor','none');
% colormap(ax2,cmocean('balance','pivot',0));
colormap(ax2,cm);

% m_gshhs_f('patch',rgb('light grey'),'edgecolor','none');
m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none');

m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
    'xtick',[],'ytick',[]);

% m_vec(0.5,-123.25,49.1088,0,0.2,...
%     'shaftwidth',2,'headlength',10,'edgeclip', 'on','headangle',30,...
%     'edgecolor','w','linewidth',1.5);
% m_text(-123.25,49.10,'0.1 m s^{-1}','fontsize',7,'horizontalalignment','center');

gr=[0.7 0.7 0.7];
m_vec(0.5,stacey.mLon,stacey.mLat,stacey.modC100U,stacey.modC100V,...
    'shaftwidth',0.5,'headlength',5,'edgeclip', 'on','headangle',30,...
    'linewidth',1,'facecolor','k','edgecolor','k');
s(1)=m_scatter(stacey.mLon,stacey.mLat,10,'filled','markerfacecolor',gr,...
    'markeredgecolor',gr);
m_vec(0.5,stacey.mLon,stacey.mLat,stacey.C100(:,1)/100,stacey.C100(:,2)/100,...
    'shaftwidth',0.5,'headlength',5,'edgeclip', 'on','headangle',30,...
    'linewidth',1,'facecolor',gr,'edgecolor',gr,'facealpha',0.5);



lnes=lines;
s(2)=m_scatter(central.lon,central.lat,40,'markerfacecolor',lnes(5,:),'markeredgecolor','none');
m_vec(0.5,central.lon,central.lat,central.modU,central.modV,...
    'shaftwidth',0.75,'headlength',10,'edgeclip', 'on','headangle',30,...
    'facecolor','k','edgecolor','k');
m_vec(0.5,central.lon,central.lat,central.mnU,central.mnV,...
    'shaftwidth',2,'headlength',7,'edgeclip', 'on','headangle',30,...
    'facecolor',lnes(5,:),'edgecolor','w','linewidth',1,'facealpha',0.5);

s(3)=m_scatter(east.lon,east.lat,40,'markerfacecolor',lnes(3,:),'markeredgecolor','none');
m_vec(0.5,east.lon,east.lat,east.modU,east.modV,...
    'shaftwidth',2,'headlength',7,'edgeclip', 'on','headangle',30,...
    'facecolor','k','edgecolor','k');
m_vec(0.5,east.lon,east.lat,east.mnU,east.mnV,...
    'shaftwidth',2,'headlength',10,'edgeclip', 'on','headangle',30,...
    'facecolor',lnes(3,:),'edgecolor','w','linewidth',1,'facealpha',0.5);

% m_annotation('textbox',[east.lon+0.05,east.lat,0,0],...
%     'string',{'East';[num2str(sqrt(east.mnU^2+east.mnV^2),'%2.2f') ' cm s^{-1}']},...
%     'fontsize',6,'horizontalalignment','center','FitBoxToText','on',...
%     'BackgroundColor','white');

s(4)=m_scatter(DDL.lon,DDL.lat,40,'markerfacecolor',lnes(4,:),'markeredgecolor','none');
m_vec(0.5,DDL.lon,DDL.lat,DDL.modU,DDL.modV,...
    'shaftwidth',2,'headlength',7,'edgeclip', 'on','headangle',30,...
    'facecolor','k','edgecolor','k');
m_vec(0.5,DDL.lon,DDL.lat,DDL.mnU,DDL.mnV,...
    'shaftwidth',2,'headlength',10,'edgeclip', 'on','headangle',30,...
    'facecolor',lnes(4,:),'edgecolor','w','linewidth',1,'facealpha',0.5);



% l=legend(s,['Central= ' num2str(sqrt(central.mnU^2+central.mnV^2),'%2.2f') ' cm s^{-1}'],...
%     ['East= ' num2str(sqrt(east.mnU^2+east.mnV^2),'%2.2f') ' cm s^{-1}'],...
%     ['DDL=' num2str(sqrt(DDL.mnU^2+DDL.mnV^2),'%2.2f') ' cm s^{-1}'],...
%     'location','northeast','fontsize',6);
text(0.925,0.1,'b)','units','normalized','fontsize',8,'Fontweight','bold');

set(findall(gcf,'-property','FontSize'),'FontSize',8)

% export_fig /ocean/sstevens/IW_project/figures/paper/vec/eulerianP1.pdf -dpdf
clearvars mean_east mean_north uv_i uv_mag ncst DDL east central

% Plot EKE map
EKE=NaN(size(mean_u));

for i=1:size(mean_u,1)
    for ii=1:size(mean_u,2)
        uprime=(squeeze(model_uv.u(i,ii,:))-mean_u(i,ii));
        vprime=(squeeze(model_uv.v(i,ii,:))-mean_v(i,ii));
        
        EKE(i,ii)=nanmean((uprime.^2+vprime.^2)/2,'all');
    end
end

MKE=(((mean_u).^2+(mean_v).^2)/2);

qlon=model_uv.lon;qlon(isnan(qlon))=0;
qlat=model_uv.lat;qlat(isnan(qlat))=0;

% Calulate percentages of region in each regime
load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;
msk=inShape(station_bound,double(qlon),double(qlat));
E2M=nansum(EKE(msk))/nansum(MKE(msk));

E2Mvec=EKE(msk)./MKE(msk);
E2Mvec(isnan(E2Mvec))=[];
msk1=E2Mvec(:)<0.75; msk2=E2Mvec(:)>=0.75 & E2Mvec(:)<=1.25; msk3=E2Mvec(:)>1.25;
ad=sum(msk1)/length(E2Mvec(:));
interm=sum(msk2)/length(E2Mvec(:));
dif=sum(msk3)/length(E2Mvec(:));

ax5=axes('Position',[0.055 0.025 0.9 0.45]);
lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -122.500002324];
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
lev=[0:0.1:1 5 10:10:90];
[CS,CH]=m_contourf(model_uv.lon,model_uv.lat,...
    log10((EKE)./(MKE)),log10(lev),'linecolor','none');
% cm=m_colmap('mBOD',length(lev)+8);cm=brighten(cm(5:end-4,:),-0.5);
cm=getPyPlot_cMap('PRGn');ct=round(length(cm)/2);
% cm=getPyPlot_cMap('seismic');ct=round(length(cm)/2);
cm(ct-18:ct+12,:)=[];

colormap(ax5,cm);
ax=colorbar('westoutside');
set(ax,'ytick',log10([0.01 0.1 1 10 100]),'yticklabel',[0 0.1 1 10 100],'tickdir','out');
ax.Position = [0.175 0.275 0.035 0.2]; % or + [left, bottom, width, height] to place it where you want

xlabel(ax,'KER','fontweight','bold');
%m_gshhs_l('patch',rgb('light grey'),'edgecolor','none');
m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none');

m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',7);

text(0.925,0.95,'c)','units','normalized','fontsize',8,'Fontweight','bold');
set(findall(gcf,'-property','FontSize'),'FontSize',8)

axes(ax2)
l=legend(s,'Stacey 1987','Central','East','DDL','location','northeast','fontsize',6);
title(l,'Node');

%%
export_fig /ocean/sstevens/IW_project/figures/paper/vec/eulerian_v7.pdf -dpdf
% export_fig /ocean/sstevens/IW_project/figures/paper/vec/eulerian_v6.jpg

%% OLD code for log color scale
% 
% load PT_AS_inlets.mat
% station_bound=alphaShape(AS.lon,AS.lat);
% station_bound.Alpha=station_bound.Alpha*2;
% msk=inShape(station_bound,double(qlon),double(qlat));
% E2M=nansum(EKE(msk))/nansum(MKE(msk));
% %%
% ax5=axes('Position',[0.64 0.15 0.4 .8150]);
% lat_lim=[50.42283746 48.252093478];
% lon_lim=[-125.82436392 -122.252345235];
% m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
% hold on
% lev=log10([0:0.25:1 10:10:10]);
% [CS,CH]=m_contourf(model_uv.lon,model_uv.lat,...
%     log10((EKE)./(MKE)),lev,'linecolor','none');
% cm=cmocean('balance','pivot',log10(1)); 
% % cm=cm(60:end-60,:);
% % tmp=round(length(cm)/2);cm(tmp-15:tmp+15,:)=[];
% colormap(ax5,cm);
% % [ax,h]=m_contfbar(ax5,[0.25 0.75],-0.0,CS,CH,'endpiece','yes','axfrac',.04,...
% %     'edgecolor','none');
% ax=colorbar;
% % set(get(ax,'ylabel'),'String','Chla (\mug/l)');
% set(ax,'ytick',log10([.5 1 2 3 5 10 100 1000]),'yticklabel',[.5 1 2 3 5 10 100 1000],...
%         'tickdir','out');
% % set(gca,'ColorScale','log')
% 
% xlabel(ax,'EKE/MKE','fontweight','bold');
% m_gshhs_l('patch',rgb('light grey'),'edgecolor','none');
% % m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none');
% m_grid('linestyle','none','linewidth',0.5,'tickdir','out','fontsize',6,...
%     'xaxisloc','top','yaxisloc','right');
% text(0.85,0.95,'c)','units','normalized','fontsize',8,'Fontweight','bold');

%% KE figs
% EKE(EKE==0)=NaN;
% MKE(MKE==0)=NaN;
% % figure('units','centimeters','outerposition',[0 0 14 14],'color','w');% lat_lim=[48.4 50.15];
% 
% ax3=axes('Position',[0.55 0.5 0.235 .47]);
% lat_lim=[50.42283746 48.252093478];
% lon_lim=[-125.82436392 -122.252345235];
% m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
% hold on
% caxis([0 5e-3]);
% [CS,CH]=m_contourf(model_uv.lon,model_uv.lat,...
%     EKE,0:15e-5:5e-3,'linecolor','none');
% colormap(ax3,mellow_rainbow);
% %m_gshhs_l('patch',rgb('light grey'),'edgecolor','none');
% m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none');
% m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
%     'xticklabel',[],'yticklabel',[]);
% text(0.8,0.05,'c)','units','normalized','fontsize',8,'Fontweight','bold');
% 
% ax4=axes('Position',[0.675 0.5 0.235 .47]);
% lat_lim=[50.42283746 48.252093478];
% lon_lim=[-125.82436392 -122.252345235];
% m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.35)
% hold on
% caxis([0 5e-3]);
% m_contourf(model_uv.lon,model_uv.lat,...
%     MKE,0:15e-5:5e-3,'linecolor','none');
% colormap(ax4,mellow_rainbow);
% [ax,h]=m_contfbar(ax4,[0.225 0.675],-0.1,CS,CH,'endpiece','yes','axfrac',.1,...
%     'edgecolor','none','levels','match','xticklabelrotation',45,'xtick',0:1e-3:4e-3);
% xlabel(ax,'KE (m^2 s^{-2})','fontweight','bold');
% 
% % m_gshhs_l('patch',rgb('light grey'),'edgecolor','none');
% m_usercoast('/ocean/rich/more/mmapbase/bcgeo/PNW.mat','patch',[1,1,1].*0.7, 'Edgecolor', 'none');
% m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
%     'xticklabel',[],'yticklabel',[]);
% text(0.8,0.05,'d)','units','normalized','fontsize',8,'Fontweight','bold');
% set(findall(gcf,'-property','FontSize'),'FontSize',8)


%% Work out transport distributions and timescale from within the Strait
load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
% station_bound=alphaShape(x,y);
station_bound.Alpha=station_bound.Alpha*2;
qlon=double(model_uv.lon);qlon(isnan(qlon))=0;
qlat=double(model_uv.lat);qlat(isnan(qlat))=0;
inside_idx=inShape(station_bound,qlon,qlat);

straitU=mean_u;straitU(~inside_idx)=NaN;
straitV=mean_v;straitV(~inside_idx)=NaN;

generalV=nanmean(straitV(straitV>0));% ms-1
generalTS=(250e3/generalV)/(60*60*24); % days

maxV=NaN(size(straitV,2),1);
mnV=NaN(size(straitV,2),1);
for i=1:size(straitV,2)
    tmpV=straitV(:,i);
    maxV(i)=max(tmpV);
    mnV(i)=nanmean(tmpV(tmpV>0)); 
    mnVp(i)=sum(~isnan(tmpV(tmpV>0)))/sum(~isnan(straitV(:,i)));
end
maxV(maxV<=0)=NaN;
maxTS=nansum(ones(length(maxV),1)*500./maxV)/(60*60*24); % days
mnTS=nansum(ones(length(mnV),1)*500./mnV)/(60*60*24); % days

%% ... and downstrait (i.e. from discovery pass.)
generalV=nanmean(straitV(straitV<0));% ms-1
generalTS=(250e3/generalV)/(60*60*24); % days

maxV=NaN(size(straitV,2),1);
mnV=NaN(size(straitV,2),1);
for i=1:size(straitV,2)
    tmpV=straitV(:,i);
    maxV(i)=min(tmpV);
    mnV(i)=nanmean(tmpV(tmpV<0)); 
end
maxV(maxV>=0)=NaN;
maxTS=nansum(ones(length(maxV),1)*500./maxV)/(60*60*24); % days
mnTS=nansum(ones(length(mnV),1)*500./mnV)/(60*60*24); % days