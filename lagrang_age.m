%% Z AGE MAPS
addpath(genpath('/ocean/sstevens'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
%%
clear

load('BCcoast');
DP=load('DP_2016age');
BP=load('BP_2016age');
lon_grid=BP.lon_grid/BP.res;lat_grid=BP.lat_grid/BP.res;
% mod_lat=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
%     'nav_lat');
% mod_lon=ncread('/results2/SalishSea/nowcast-green.201812/01apr15/SalishSea_1d_20150401_20150401_grid_T.nc',...
%     'nav_lon');

mod_lat=ncread('/ocean/sallen/allen/research/MEOPAR/grid/mesh_mask201702.nc',...
'nav_lat');
mod_lon=ncread('/ocean/sallen/allen/research/MEOPAR/grid/mesh_mask201702.nc',...
'nav_lon');

fname='/ocean/rich/more/mmapbase/noaa_bc3/british_columbia_3_msl_2013.nc';

lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
lat_lim=[48.09865765 50.211000657653];
lon_lim=[-125.356524173 -122.2500002324];

ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
lnes=lines;
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);
[Xz,Yz]=meshgrid(lat(ilat),lon(ilon));

FS=load('AS_Fraser.mat');
station_bound=alphaShape(FS.X,FS.Y);
msk=~inShape(station_bound,Yz,Xz);
Z(msk)=NaN;
%%
DP.count_map=NaN(size(lat_grid));
DP.mean_map=NaN(size(lat_grid));
DP.endpoint_counter=NaN([size(lat_grid),9]);
DP.std_map=NaN(size(lat_grid));
DP.zmean_map=NaN(size(DP.zage_map));
DP.zendpoint_counter=NaN([size(DP.zage_map),size(DP.zage_map,3)]);

BP.count_map=NaN(size(lat_grid));
BP.mean_map=NaN(size(lat_grid));
BP.endpoint_counter=NaN([size(lat_grid),9]);
BP.std_map=NaN(size(lat_grid));
BP.zmean_map=NaN(size(DP.zage_map));
BP.zendpoint_counter=NaN([size(DP.zage_map),size(DP.zage_map,3)]);

for i=1:size(DP.age_map,1)
    for ii=1:size(DP.age_map,2)
        DP.count_map(i,ii)=length(DP.age_map{i,ii});
        DP.mean_map(i,ii)=mean(DP.age_map{i,ii});
        DP.std_map(i,ii)=std(DP.age_map{i,ii});
        
        BP.count_map(i,ii)=length(BP.age_map{i,ii});
        BP.mean_map(i,ii)=mean(BP.age_map{i,ii});
        BP.std_map(i,ii)=std(BP.age_map{i,ii});
        
        for iii=1:size(DP.endpoint_map,3)
            DP.endpoint_counter(i,ii,iii)=length(DP.endpoint_map{i,ii,iii});
            BP.endpoint_counter(i,ii,iii)=length(BP.endpoint_map{i,ii,iii});
        end
      
    end
end


for i=1:size(DP.zage_map,1)
    for ii=1:size(DP.zage_map,2)
        
        if length(DP.zage_map{i,ii})>2
            DP.zmean_map(i,ii)=mean(DP.zage_map{i,ii});
        end
        
        if length(BP.zage_map{i,ii})>20
            BP.zmean_map(i,ii)=mean(BP.zage_map{i,ii});
        end

        for iii=1:size(DP.zendpoint_map,3)
            DP.zendpoint_counter(i,ii,iii)=length(DP.zendpoint_map{i,ii,iii});
            BP.zendpoint_counter(i,ii,iii)=length(BP.zendpoint_map{i,ii,iii});
        end
    end
end

idx_south=lon_grid<-124.5 & lat_grid<49;
idx_north=lon_grid<-125.5 & lat_grid>50;

blocklon=BlockMean(DP.lon_grid/100,4,4);
blocklat=BlockMean(DP.lat_grid/100,4,4);

% Smoothing
filt1=ones(4,4);
for i=1:size(DP.endpoint_counter,3)
    
    tmp_dp=squeeze(DP.endpoint_counter(:,:,i));
    DP.left_north(i)=nansum(tmp_dp(idx_north));
    DP.left_south(i)=nansum(tmp_dp(idx_south));
    tmp_dp(logical(idx_north+idx_south))=0;
    DP.remain(i)=nansum(tmp_dp(:));
    
    tmp_bp=squeeze(BP.endpoint_counter(:,:,i));
    BP.left_north(i)=nansum(tmp_bp(idx_north));
    BP.left_south(i)=nansum(tmp_bp(idx_south));
    tmp_bp(logical(idx_north+idx_south))=0;
    BP.remain(i)=nansum(tmp_bp(:));

    DP.endpoint_smooth(:,:,i)=conv2(tmp_dp,filt1,'same');
    BP.endpoint_smooth(:,:,i)=conv2(tmp_bp,filt1,'same');
    
    % normalise endpoint map BY THE ?? PERCENTILE OF THE DATA
    DP.endpoint_smooth(:,:,i)=DP.endpoint_smooth(:,:,i)./...
        prctile(DP.endpoint_smooth(:,:,i),99.9,'all');
    BP.endpoint_smooth(:,:,i)=BP.endpoint_smooth(:,:,i)./...
        prctile(BP.endpoint_smooth(:,:,i),99.9,'all');
    
    % Try block mean
    DP.endpoint_block(:,:,i)=BlockMean(tmp_dp,4,4);
    BP.endpoint_block(:,:,i)=BlockMean(tmp_bp,4,4);
        
    DP.endpoint_block(:,:,i)=DP.endpoint_block(:,:,i)./...
        prctile(DP.endpoint_block(:,:,i),100,'all');
    BP.endpoint_block(:,:,i)=BP.endpoint_block(:,:,i)./...
        prctile(BP.endpoint_block(:,:,i),100,'all');
end

[X,Y]=meshgrid(1:450,lat_grid(:,1));X=rot90(X,3);Y=fliplr(rot90(Y,3));
load PT_bottom
idx=Y(1,:)<min(a) | Y(1,:)>max(a);bottomY=Y(1,~idx);
bottom=interp1(a,b,bottomY,'spline');

% Smoothing
filt1=ones(5,5);

blockX=BlockMean(X,5,2);
blockY=BlockMean(Y,5,2);

for i=1:size(DP.zendpoint_counter,3)
    tmp_dp=squeeze(DP.zendpoint_counter(:,:,i));
    tmp_bp=squeeze(BP.zendpoint_counter(:,:,i));
    
    DP.zendpoint_smooth(:,:,i)=conv2(tmp_dp,filt1,'same');
    BP.zendpoint_smooth(:,:,i)=conv2(tmp_bp,filt1,'same');
    
    % Try block mean
    DP.zendpoint_block(:,:,i)=BlockMean(tmp_dp,5,2);
    BP.zendpoint_block(:,:,i)=BlockMean(tmp_bp,5,2);
    
    % normalise endpoint map BY THE ?? PERCENTILE OF THE DATA
    DP.zendpoint_smooth(:,:,i)=DP.zendpoint_smooth(:,:,i)./...
        prctile(DP.zendpoint_smooth(:,:,i),100,'all');
    BP.zendpoint_smooth(:,:,i)=BP.zendpoint_smooth(:,:,i)./...
        prctile(BP.zendpoint_block(:,:,i),100,'all');
    
    DP.zendpoint_block(:,:,i)=DP.zendpoint_block(:,:,i)./...
        prctile(DP.zendpoint_block(:,:,i),100,'all');
    BP.zendpoint_block(:,:,i)=BP.zendpoint_block(:,:,i)./...
        prctile(BP.zendpoint_block(:,:,i),100,'all');
    
end

DP.endpoint_smooth(DP.endpoint_smooth>1)=1;
BP.endpoint_smooth(BP.endpoint_smooth>1)=1;
DP.zendpoint_smooth(DP.zendpoint_smooth>1)=1;
BP.zendpoint_smooth(BP.zendpoint_smooth>1)=1;
DP.endpoint_block(DP.endpoint_block>1)=1;
BP.endpoint_block(BP.endpoint_block>1)=1;
DP.zendpoint_block(DP.zendpoint_block>1)=1;
BP.zendpoint_block(BP.zendpoint_block>1)=1;
%%
% Remove mean age from release region to better represent transit times
idx=dsearchn([BP.lon_grid(:)/100 BP.lat_grid(:)/100],...
    [-122.9962 48.8496]);
[row,col]=ind2sub(size(BP.lon_grid),idx);

BP.zmean_map=BP.zmean_map-...
    nanmean(BP.mean_map(row-2:row+2,col-2:col+2),'all');
BP.mean_map=BP.mean_map-...
    nanmean(BP.mean_map(row-2:row+2,col-2:col+2),'all');

idx=dsearchn([DP.lon_grid(:)/100 DP.lat_grid(:)/100],...
    [mod_lon(124,737) mod_lat(124,737)]);
[row,col]=ind2sub(size(DP.lon_grid),idx);

DP.zmean_map=DP.zmean_map-...
    nanmean(DP.mean_map(row-2:row+2,col-2:col+2),'all');
DP.mean_map=DP.mean_map-...
    nanmean(DP.mean_map(row-2:row+2,col-2:col+2),'all');

% save /ocean/sstevens/IW_project/data/lagrangianM.mat DP BP

%% Plot
% Boundary Pass
figure('units','centimeters','outerposition',[0 0 19 22],'color','w');
ax1=axes('Position',[0.055 0.45 0.4 0.45]);

m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');
[CS,CH]=m_contourf(lon_grid,lat_grid,BP.mean_map,-20:20:260,'linestyle','none');
colormap(ax1,mod_blue_orange);
caxis([0 160]);
[ax,h]=m_contfbar(ax1,[0.15 0.95],1.07,CS,0:20:160,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',0:20:160,'XTickLabelRotation',30,'levels','match');
ylabel(ax,{'Transit time';'(days)'},'fontsize',8,'color','k','fontweight','bold','rotation',0);

[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,RGB,'edgecolor','none');
end
text(1.025,0.05,'a)','fontsize',8,'Fontweight','bold','units','normalized');

m_scatter(mod_lon(305,350),mod_lat(305,350),40,'filled','w');
m_scatter(mod_lon(305,350),mod_lat(305,350),30,'x','k','linewidth',1);

m_line([-124.4331;-123.6],[49.1633;49.7210],'color','w','linestyle','--');
m_line([-123.6413;-122.8],[48.5701;49.0969],'color','w','linestyle','--');

m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',8);

axpos={[0.5 0.45 0.1 0.45];[0.62 0.45 0.1 0.45];...
    [0.74 0.45 0.1 0.45];[0.86 0.45 0.1 0.45]};

lats_lim=[50.32283746 48.52093478];
lons_lim=[-125.62436392 -122.252345235];
cnum=[1 2 4 8];
labstr={'b)','c)','d)','e)'};
for i=1:4
    ax(i)=axes('Position',axpos{i});
    %     ax(i)=subplot(4,2,4+i);
    m_proj('oblique','lon',lons_lim,'lat',lats_lim,'dir','vert','aspect',0.225)
    hold on
    colormap(ax(i),cmocean('amp'));
    [CS,CH]=m_contourf(blocklon,blocklat,BP.endpoint_block(:,:,cnum(i)),...
        [0.01 0.1:0.1:1],'linestyle','none');
    caxis([0.01 1]);
    
    m_gshhs_f('patch',rgb('light grey'),'edgecolor','none');
    m_contour(Yz,Xz,Z,[-10 -10],'linecolor',rgb('light grey')-0.2,'linewidth',0.5);

    if i==1
        m_northarrow(-125.5273,50.0613,0.2,'type',2,'linewidth',0.5);
        m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
            'xaxisloc','bottom','yaxisloc','right','fontsize',6,'ytick',[49 50]); 
    else
        m_grid('linestyle','none','tickdir','out','xticklabel',[],'yticklabels',[]);
    end
    
    m_scatter(mod_lon(305,350),mod_lat(305,350),70,'x','k','linewidth',3);
    text(0.75,0.05,labstr{i},'fontsize',8,'Fontweight','bold','units','normalized');
    title(['Day ' num2str(cnum(i)*20)],'fontsize',8,'Fontweight','normal');
end

[ax(i+1),h]=m_contfbar(ax(1),[0.95 4.5],1.08,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'XTickLabelRotation',30,'xtick',[0.01 0.1:0.1:1]);
ylabel(ax(i+1),{'Particle';'Concentration'},'fontsize',8,'color','k','fontweight','bold','rotation',0);

%%%%%%%%%%%%%%%%%%%%
% Depth-Latitude plot

[X,Y]=meshgrid(1:450,lat_grid(:,1));X=rot90(X,3);Y=fliplr(rot90(Y,3));

ax1=axes('position',[0.075,0.15,0.5,0.25],'fontsize',8);
set(gca,'color','k');
hold on
RGB=rgb('light grey');
[CS,CH]=contourf(Y,X,BP.zmean_map,-20:20:160,'linestyle','none');
caxis([0 160]);

tlab=48.9:0.2:50;
for i=1:length(tlab)
    tlabt{i}=[num2str(tlab(i)) char(176)];
end
set(gca,'ydir','reverse','tickdir','out','xtick',tlab,...
    'xticklabels',tlabt);
colormap(ax1,mod_blue_orange);
xlim([48.74 50.01]);

area(bottomY,bottom,450,'FaceColor',rgb('light grey'),'edgecolor','none');
plot(ones(2,1)*mean([49.1633 49.7210]),[0 interp1(bottomY,bottom,...
    mean([49.1633 49.7210]))],'w--','linewidth',2);
plot(ones(2,1)*mean([48.5701;49.0969]),[0 interp1(bottomY,bottom,...
    mean([48.5701;49.0969]))],'w--','linewidth',2);
ylim([0 450]); xlim([min(bottomY) max(bottomY)]);
text(49.9,400,'f)','fontsize',8,'Fontweight','bold');
text(mean([49.1633 49.7210]),-15,'NB','fontsize',8,...
    'HorizontalAlignment','center');
text(mean([48.5701;49.0969]),-15,'SB','fontsize',8,...
    'HorizontalAlignment','center');
scatter([48.8176 48.8176 48.8176],[50 75 100],40,'o','w','filled',...
    'markerfacealpha',1,'markeredgecolor','none');
scatter([48.8176 48.8176 48.8176],[50 75 100],30,'x','k','linewidth',2);

ylabel('Depth (m)','fontsize',8,'Fontweight','bold');
xlabel('Latitude (^\circN)','Fontweight','bold','fontsize',8);

cnum=[1 2 4 8];% 12 size(DP.endpoint_smooth,3)];
labstr={'g)','h)','i)','j)','k)','l)'};
psp=[[linspace(0.625,0.81,2)';linspace(0.625,0.81,2)'] ...
    [ones(2,1)*0.35;ones(2,1)*0.2] ...
    ones(4,1)*0.15 ones(4,1)*0.1];
psp(:,2)=psp(:,2)-0.05;

for i=1:length(cnum)
    ax(i)=axes('Position',psp(i,:),'fontsize',8);
    hold on
    [CS,CH]=contourf(blockY,blockX,BP.zendpoint_block(:,:,cnum(i)),...
        [0.01 0.1:0.1:1],'linestyle','none');
    colormap(ax(i),cmocean('amp'));
    if i<=2
        set(gca,'ydir','reverse','tickdir','out','xtick',[]);
    else
        set(gca,'ydir','reverse','tickdir','out','xtick',tlab(1:2:end),...
            'xticklabels',tlabt(1:2:end),'xticklabelrotation',30);
    end
    
    if rem(i,2)==0
        set(gca,'ytick',[]);
    end
        
    ylim([0 450]); xlim([min(bottomY) max(bottomY)]);
    area(bottomY,bottom,450,'FaceColor',rgb('light grey'),'edgecolor','none');
    plot(ones(2,1)*mean([49.1633 49.7210]),[0 interp1(bottomY,bottom,...
        mean([49.1633 49.7210]))],'--','linewidth',1,'color',[1 1 1]*0.0);
    plot(ones(2,1)*mean([48.5701;49.0969]),[0 interp1(bottomY,bottom,...
        mean([48.5701;49.0969]))],'--','linewidth',1,'color',[1 1 1]*0.0);
    title(['Day ' num2str(cnum(i)*20)],'fontsize',8,'fontweight','normal');
    text(49.8,390,labstr{i},'fontsize',8,'Fontweight','bold');
    scatter([48.8176 48.8176 48.8176],[50 75 100],5,'x','k','linewidth',1);
end

%%
%export_fig '/ocean/sstevens/IW_project/figures/paper/zage_map_BP2.png' -png -m3
export_fig '/ocean/sstevens/IW_project/figures/paper/vec/zage_map_BP2.pdf' -dpdf
%%
% Discovery Pass
figure('units','centimeters','outerposition',[0 0 19 22],'color','w');
ax1=axes('Position',[0.055 0.45 0.4 0.45]);

m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');
[CS,CH]=m_contourf(lon_grid,lat_grid,DP.mean_map,-20:20:160,'linestyle','none');
colormap(ax1,mod_blue_orange);
caxis([0 160]);
[ax,h]=m_contfbar(ax1,[0.15 0.95],1.07,CS,0:20:160,'endpiece','no','axfrac',.03,'fontsize',8,...
    'xtick',0:20:160,'XTickLabelRotation',30,'levels','match');
ylabel(ax,{'Mean Age';'(days)'},'fontsize',8,'color','k','fontweight','bold','rotation',0);

[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,RGB,'edgecolor','none');
end
text(1.025,0.05,'a)','fontsize',8,'Fontweight','bold','units','normalized');

m_scatter(mod_lon(124,737),mod_lat(124,737),40,'filled','w');
m_scatter(mod_lon(124,737),mod_lat(124,737),30,'x','k','linewidth',1);

m_line([-124.4331;-123.6],[49.1633;49.7210],'color','w','linestyle','--');
m_line([-123.6413;-122.8],[48.5701;49.0969],'color','w','linestyle','--');

m_grid('linestyle','none','linewidth',0.5,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','left','fontsize',8);

axpos={[0.5 0.45 0.1 0.45];[0.62 0.45 0.1 0.45];...
    [0.74 0.45 0.1 0.45];[0.86 0.45 0.1 0.45]};

lats_lim=[50.32283746 48.52093478];
lons_lim=[-125.62436392 -122.252345235];
cnum=[1 2 4 8];
labstr={'b)','c)','d)','e)'};
for i=1:4
    ax(i)=axes('Position',axpos{i});
    %     ax(i)=subplot(4,2,4+i);
    m_proj('oblique','lon',lons_lim,'lat',lats_lim,'dir','vert','aspect',0.225)
    hold on
    colormap(ax(i),cmocean('amp'));
    [CS,CH]=m_contourf(blocklon,blocklat,DP.endpoint_block(:,:,cnum(i)),...
        [0.01 0.1:0.1:1],'linestyle','none');
    caxis([0.01 1]);
    
    m_gshhs_f('patch',rgb('light grey'),'edgecolor','none');
    m_contour(Yz,Xz,Z,[-10 -10],'linecolor',rgb('light grey')-0.2,'linewidth',0.5);

    if i==1
        m_northarrow(-125.5273,50.0613,0.2,'type',2,'linewidth',0.5);
        m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
            'xaxisloc','bottom','yaxisloc','right','fontsize',6,'ytick',[49 50]); 
    else
        m_grid('linestyle','none','tickdir','out','xticklabel',[],'yticklabels',[]);
    end
    
    m_scatter(mod_lon(124,737),mod_lat(124,737),70,'x','k','linewidth',3);
    text(0.75,0.05,labstr{i},'fontsize',8,'Fontweight','bold','units','normalized');
    title(['Day ' num2str(cnum(i)*20)],'fontsize',8,'Fontweight','normal');
end

[ax(i+1),h]=m_contfbar(ax(1),[0.95 4.5],1.08,CS,CH,'endpiece','no','axfrac',.03,'fontsize',8,...
    'XTickLabelRotation',30,'xtick',[0.01 0.1:0.1:1]);
ylabel(ax(i+1),{'Particle';'Concentration'},'fontsize',8,'color','k','fontweight','bold','rotation',0);

%%%%%%%%%%%%%%%%%%%%
% Depth-Latitude plot

[X,Y]=meshgrid(1:450,lat_grid(:,1));X=rot90(X,3);Y=fliplr(rot90(Y,3));

ax1=axes('position',[0.075,0.15,0.5,0.25],'fontsize',8);
set(gca,'color','k');
hold on
RGB=rgb('light grey');
[CS,CH]=contourf(Y,X,DP.zmean_map,-20:20:160,'linestyle','none');
caxis([0 160]);

tlab=48.9:0.2:50;
for i=1:length(tlab)
    tlabt{i}=[num2str(tlab(i)) char(176)];
end
set(gca,'ydir','reverse','tickdir','out','xtick',tlab,...
    'xticklabels',tlabt);
colormap(ax1,mod_blue_orange);
xlim([48.74 50.01]);

area(bottomY,bottom,450,'FaceColor',rgb('light grey'),'edgecolor','none');
plot(ones(2,1)*mean([49.1633 49.7210]),[0 interp1(bottomY,bottom,...
    mean([49.1633 49.7210]))],'w--','linewidth',2);
plot(ones(2,1)*mean([48.5701;49.0969]),[0 interp1(bottomY,bottom,...
    mean([48.5701;49.0969]))],'w--','linewidth',2);
ylim([0 450]); xlim([min(bottomY) max(bottomY)]);
text(49.9,400,'f)','fontsize',8,'Fontweight','bold');
text(mean([49.1633 49.7210]),-15,'NB','fontsize',8,...
    'HorizontalAlignment','center');
text(mean([48.5701;49.0969]),-15,'SB','fontsize',8,...
    'HorizontalAlignment','center');
scatter(49.95,50,40,'o','w','filled',...
    'markerfacealpha',1,'markeredgecolor','none');
scatter(49.95,50,30,'x','k','linewidth',2);

ylabel('Depth (m)','fontsize',8,'Fontweight','bold');
xlabel('Latitude (^\circN)','Fontweight','bold','fontsize',8);

cnum=[1 2 4 8];% 12 size(DP.endpoint_smooth,3)];
labstr={'g)','h)','i)','j)','k)','l)'};
psp=[[linspace(0.625,0.81,2)';linspace(0.625,0.81,2)'] ...
    [ones(2,1)*0.35;ones(2,1)*0.2] ...
    ones(4,1)*0.15 ones(4,1)*0.1];
psp(:,2)=psp(:,2)-0.05;

for i=1:length(cnum)
    ax(i)=axes('Position',psp(i,:),'fontsize',8);
    hold on
    [CS,CH]=contourf(blockY,blockX,DP.zendpoint_block(:,:,cnum(i)),...
        [0.01 0.1:0.1:1],'linestyle','none');
    colormap(ax(i),cmocean('amp'));
    if i<=2
        set(gca,'ydir','reverse','tickdir','out','xtick',[]);
    else
        set(gca,'ydir','reverse','tickdir','out','xtick',tlab(1:2:end),...
            'xticklabels',tlabt(1:2:end),'xticklabelrotation',30);
    end
    
    if rem(i,2)==0
        set(gca,'ytick',[]);
    end
        
    ylim([0 450]); xlim([min(bottomY) max(bottomY)]);
    area(bottomY,bottom,450,'FaceColor',rgb('light grey'),'edgecolor','none');
    plot(ones(2,1)*mean([49.1633 49.7210]),[0 interp1(bottomY,bottom,...
        mean([49.1633 49.7210]))],'--','linewidth',1,'color',[1 1 1]*0.0);
    plot(ones(2,1)*mean([48.5701;49.0969]),[0 interp1(bottomY,bottom,...
        mean([48.5701;49.0969]))],'--','linewidth',1,'color',[1 1 1]*0.0);
    title(['Day ' num2str(cnum(i)*20)],'fontsize',8,'fontweight','normal');
    text(49.8,390,labstr{i},'fontsize',8,'Fontweight','bold');
    scatter(49.95,50,30,'x','k','linewidth',1);
end


%%
% export_fig '/ocean/sstevens/IW_project/figures/paper/zage_map_DP2.png' -png -m3
% export_fig '/ocean/sstevens/IW_project/figures/paper/vec/zage_map_DP2.pdf' -dpdf
