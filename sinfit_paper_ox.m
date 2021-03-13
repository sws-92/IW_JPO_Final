%% Fit sinusoids to all stations and use objective mapping to to create map
% of coldest yearday for entire strait

% Stations to use
% Stratogem
% PSF
% Hakai
% IOS
% any Nanoose?

clear
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));

%%
clear
load /ocean/sstevens/IW_project/data/all_ctd.mat
load /ocean/sstevens/IW_project/data/no_inlet_msk.mat
load /ocean/sstevens/IW_project/data/thalweg.mat
lnes=lines;
grey=rgb('light grey');

% % Remove IoS stations
% clc
% q=input('Would you like to include IOS stations in kriging? y/n: ','s');
% if strcmp(q,'n')
%     idx=strcmp(all_ctd.dataset,'IOS_ctd');
%     all_ctd.lat(idx)=[];
%     all_ctd.lon(idx)=[];
%     all_ctd.ox(:,idx)=[];
%     all_ctd.station(idx)=[];
%     all_ctd.dataset(idx)=[];
% end

% make histogram of all stations
[sort_stns,sort_idx]=sort(all_ctd.station);
dataset_ref=all_ctd.dataset(sort_idx);
[uniq_stns,uniq_idx]=unique(sort_stns);
dataset_ref=dataset_ref(uniq_idx);
[uniq_idx,ref_idx]=sort(uniq_idx);
dataset_ref=dataset_ref(ref_idx);

for i=1:length(uniq_stns)
    if i~=length(uniq_stns)
        num_stns(uniq_idx(i):uniq_idx(i+1)-1)=i;
    else
        num_stns(uniq_idx(i):length(sort_stns))=i;
    end
end

[hist_station,edges]=histcounts(num_stns,1:length(uniq_stns)+1);

% Need to remove stations with less than 9 months of samples
hidx=hist_station>10;
rel_stations=hist_station(hidx);
rel_dataset=dataset_ref(hidx);
rel_station_names=uniq_stns(hidx);
[rel_dataset,didx]=sort(rel_dataset);
rel_stations=rel_stations(didx);
rel_station_names=rel_station_names(didx);
celld=unique(rel_dataset);


%% Create sinusoid fits to find coldest yearday

% Loop through all stations with more that 10 profiles
min_day=NaN(1,length(rel_station_names));
fitted_harms=NaN(365,length(rel_station_names));
station_lat=NaN(1,length(length(rel_station_names)));
station_lon=NaN(1,length(length(rel_station_names)));
bad_switch=NaN(1,length(length(rel_station_names)));
oamp=NaN(1,length(rel_station_names));
oxmean=NaN(1,length(rel_station_names));
h=NaN(1,length(length(rel_station_names)));
pval=NaN(1,length(length(rel_station_names)));
min_dayCI=NaN(length(rel_station_names),2);
oampCI=NaN(length(rel_station_names),2);
min_daySD=NaN(length(rel_station_names),1);
oampSD=NaN(length(rel_station_names),1);

q='n';
% q=input('Do you want to plot stations and sinusoids? y/n: ','s');
% if strcmp(q,'y')
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     ax3=axes('Position',[.05 0.05 .9 .7]);
%     thalweg_map(0,0,0,0);
%     ax2=axes('Position',[.15 .75 .5 .2]);
% end

for i=1:length(rel_station_names)
    
    % Store lat & lon
    idx=strcmp(all_ctd.station,rel_station_names{i});
    station_lon(i)=nanmean(all_ctd.lon(idx));
    station_lat(i)=nanmean(all_ctd.lat(idx));
    
    % Pull relevant station data from structure and find yearday
    stn_idx=strcmp(all_ctd.station,rel_station_names{i});
    tmp_mtime=all_ctd.mtime(stn_idx);
    nan_msk=isnan(tmp_mtime);
    tmp_mtime=tmp_mtime(~nan_msk);
    tmp_day=day(datetime(datestr(tmp_mtime)),'dayofyear');
    [tmp_day,sidx]=sort(tmp_day);
    tmp_ox=nanmean(all_ctd.ox(90:110,stn_idx));
    tmp_ox(nan_msk)=[];
    tmp_ox=tmp_ox(sidx);
    tmp_ox=inpaint_nans(tmp_ox');
    
    % fit sinusoid
    
    [oamp(i),phase,frac,offset,yy]=fit_harmonics(tmp_ox',tmp_day',...
        1,365,0.01);
    yt=offset;
    for tt=1:365
        fitted_harm(tt)=offset+oamp(i)*cos(2*pi*tt/365 + phase);
    end
    
    [~,min_day(i)]=min(fitted_harm);
    fitted_harms(:,i)=fitted_harm;
    oxmean(i)=nanmean(fitted_harm);
    
    % Check goodness-of-fit
    h(i)=nanmean(abs(tmp_ox-interp1(1:365,fitted_harm,tmp_day)));
    
    % Bootstrapping
    % Create sample indices
    I=bootrnd(length(tmp_ox),500);
    BSoamp=NaN(length(tmp_ox),1);
    BSphase=NaN(length(tmp_ox),1);
    BSmin_day=NaN(length(tmp_ox),1);
    
    % Loop through bootstrap samples and do harmonic analysis
    for ii=1:size(I,2)
        sample=sortrows([tmp_ox(I(:,ii)) tmp_day(I(:,ii))],2);
        [BSoamp(ii),BSphase]=fit_harmonics(sample(:,1),sample(:,2),...
            1,365,0.01);
        for tt=1:365
            BSfitted_harm(tt)=offset+BSoamp(ii)*cos(2*pi*tt/365 + BSphase);
        end
        [~,BSmin_day(ii)]=min(BSfitted_harm);
    end
    
    % Find bootstrap 95% CI
    SEM=std(BSoamp)/sqrt(length(BSoamp));
    ts = tinv([0.025  0.975],length(BSoamp)-1);
    oampCI(i,:)=nanmean(BSoamp)+ts*SEM;
    oampSD(i)=nanstd(BSoamp);
    
    SEM=std(BSmin_day)/sqrt(length(BSmin_day));
    ts = tinv([0.025  0.975],length(BSmin_day)-1);
    min_dayCI(i,:)=nanmean(BSmin_day)+ts*SEM;
    min_daySD(i)=nanstd(BSmin_day);
%     
%     
%     if oampSD(i)>30 %   strcmp(q,'y') && min_day(i)~=1 || strcmp(rel_station_names(i),'IS-1')
%         axes(ax3);
%         s=m_scatter(station_lon(i),station_lat(i),50,'filled');
%         axes(ax2)
%         p(1)=plot(1:365,fitted_harm);
%         hold on
%         p(2)=scatter(tmp_day,tmp_ox);
%         p(3)=scatter(min_day(i),fitted_harm(min_day(i)),[],'b','filled');
%         l=legend(rel_dataset{i});
%         disp(rel_station_names{i});
%         q2=input('Is this a good station? y/n: ','s');
%         if strcmp(q2,'n')
%             bad_switch(idx)=1;
%         end
%         delete(s)
%         delete(p)
%         delete(l)
%     end
%     
    % Save certain stations for plotting elsewhere
    if strcmp(rel_station_names(i),'PR-2')
        PR2.oamp=oamp(i); PR2.fitted_harm=fitted_harm;
        PR2.day=tmp_day; PR2.ox=tmp_ox;
        PR2.lat=station_lat(i);PR2.lon=station_lon(i);
        PR2.min_day=min_day(i);
    elseif strcmp(rel_station_names(i),'GO-4')
        GO4.oamp=oamp(i); GO4.fitted_harm=fitted_harm;
        GO4.day=tmp_day; GO4.ox=tmp_ox;
        GO4.lat=station_lat(i);GO4.lon=station_lon(i);
        GO4.min_day=min_day(i);
    end

end

% Remove broken stations
idx=min_day==1 | min_day<100 | oamp>145 | oampSD'>30;

% Remove Texada Island, Saanich, Comox stations, and bad BS stations
idx_TI_Saan=contains(rel_station_names,{'IOS_22';'LS-3';'LS-4';'LS-5';'BS-6';...
    'IS-1'});
idx=logical(idx+idx_TI_Saan);

min_day(idx)=[]; station_lat(idx)=[]; station_lon(idx)=[];
rel_dataset(idx)=[]; rel_station_names(idx)=[]; fitted_harms(:,idx)=[];
rel_stations(idx)=[]; oamp(idx)=[];oxmean(idx)=[]; oampSD(idx)=[];
oampCI(idx,:)=[]; min_daySD(idx)=[]; min_dayCI(idx,:)=[];

% Save data for plotting elsewhere
save('/ocean/sstevens/IW_project/data/represent_stns_ox.mat','PR2','GO4');

[slope,b]=TheilSen([min_day' oxmean']);

%% set up kriging of coldest yearday
load /ocean/sstevens/IW_project/data/thalweg.mat
load('BCcoast');

% create grids and lines for kriging
west_line=linspace(-125.9,-123.5,1000);
east_line=linspace(-124.4,-122,1000);
center_line=[linspace(-125.15,-122.75,1000);linspace(50.15,48.4,1000)];

% center_line_dist=NaN(1,1000);
% for i=1:length(center_line)
%     center_line_dist(i)=gsw_distance([center_line(1,i) center_line(1,end)],...
%         [center_line(2,i) center_line(2,end)]);
%     i
% end
load centre_line_dist.mat

grid_lat=repmat([linspace(50.15,48.4,1000)]',1,1000);
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

lat_lim=[48.4 50.15];
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

% Find station distances along center line
distance_idx=NaN(1,length(station_lat));
for i=1:length(station_lat)
    [~,idxb]=min(abs(center_line(2,:)-station_lat(i)));
    distance_N(i)=center_line_dist(idxb);
end

% Sort by distance
[distance_N,Nidx]=sort(distance_N);
Dsort_min_day=min_day(Nidx);
Dsort_oxmean=oxmean(Nidx);
Dsort_oamp=oamp(Nidx);
Dsort_harms=fitted_harms(:,Nidx);
Dsort_lat=station_lat(Nidx);
Dsort_lon=station_lon(Nidx);
Dsort_oampSD=oampSD(Nidx);
Dsort_min_daySD=min_daySD(Nidx);

save('/ocean/sstevens/IW_project/data/sinfit_ox.mat');

%% Plotting

% Set up XY projection
lat_limV=[min(grid_lat(:)) max(grid_lat(:))];
lon_limV=[min(grid_lon(:)) max(grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(station_lon,station_lat);
[Xgrid,Ygrid]=m_ll2xy(grid_lon,grid_lat);

dday = variogram([X' Y'],min_day','maxdist',49000,'plotit',false);

[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',min_day',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25],...
    [station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -122.500002324];

celld_label=strrep(celld,'_ctd','');
celld_label=strrep(celld_label,'Noos','Nanoose');

figure('units','centimeters','outerposition',[0 0 17 25]);
ax3=axes('Position',[0.05 0.55 1 0.45]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(grid_lon,grid_lat,zi-220,[0:10:100],'linecolor','k');
colormap(ax3,m_colmap('jet'));
m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');

[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,RGB);
end

% markers='^dosh';
% for i=1:length(celld)
%     [X,Y]=m_ll2xy(station_lon(strcmp(rel_dataset,celld{i})),...
%         station_lat(strcmp(rel_dataset,celld{i})));
%     s(i)=scatter(X,Y,10,markers(i),'k','filled');
% end

s=m_scatter(station_lon,station_lat,5,'x','k');
m_text(-122.75,50,'a)','fontweight','bold','fontsize',8);

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[ax,h]=m_contfbar(ax3,0.11,[0.0 0.95],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
ylabel(ax,{'Age (days)'},'fontsize',8,'color','k','fontweight','bold');

% Plot semivariogram
axes('Position',[0.255 0.56 0.28 0.2]);
set(gca, 'box','on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
text(0.1,0.1,'b)','units','normalized','fontweight','bold','fontsize',8);
axes('Position',[0.32 0.595 0.2 0.15]);
variogramfit_plot(dday.distance/1000,dday.val,[],[],[],'plotit',true);
set(gca,'fontsize',8);

% kriging of amplitude %
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(station_lon,station_lat);

damp = variogram([X' Y'],oamp',...
    'maxdist',95000,'plotit',false);

[a,c,n,vstruct] = variogramfit(damp.distance,damp.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',oamp',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25;],...
    [station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

% Set up figure
lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -122.500002324];

% Set up figure
ax1=axes('Position',[0.05 0.05 1 0.45]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(grid_lon,grid_lat,zi,[10:5:30 40:10:70]);
colormap(ax1,flipud(m_colmap('jet')));
m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');

[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,RGB);
end

s=m_scatter(station_lon,station_lat,5,'x','k');
m_text(-122.75,50,'c)','fontweight','bold','fontsize',8);

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[ax,h]=m_contfbar(ax1,0.11,[0.0 0.95],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
ylabel(ax,{'Amplitude (\mumol kg^{-1})'},'fontsize',8,'color','k','fontweight','bold');

% Plot semivariogram
axes('Position',[0.255 0.06 0.28 0.2]);
set(gca, 'box','on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
text(0.1,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);
axes('Position',[0.32 0.095 0.2 0.15]);
variogramfit_plot(damp.distance/1000,damp.val,[],[],[],'plotit',true);
set(gca,'fontsize',8);


%%
set(gcf, 'Color', 'w');
export_fig /ocean/sstevens/IW_project/figures/paper/ox_minday_amp.png -png -m3

%% Plot sinusoids

% Find station distances along center line
distance_idx=NaN(1,length(station_lat));
for i=1:length(station_lat)
    [~,idxb]=min(abs(center_line(2,:)-station_lat(i)));
    distance_N(i)=center_line_dist(idxb);
end

% Sort by distance
[distance_N,Nidx]=sort(distance_N);
Dsort_min_day=min_day(Nidx);
Dsort_oxmean=oxmean(Nidx);
Dsort_oamp=oamp(Nidx);
Dsort_harms=fitted_harms(:,Nidx);
Dsort_lat=station_lat(Nidx);
Dsort_lon=station_lon(Nidx);
Dsort_oampSD=oampSD(Nidx);
Dsort_min_daySD=min_daySD(Nidx);

% Longest timescale? Southernmost station to maximum (which is NW of
% Texada)
templong=max(min_day)-Dsort_min_day(1);
disp(['Inflow to NW time (minday)=',num2str(templong)]);

% create colmap
mycol=m_colmap('jet',length(Dsort_min_day));

% Set up figure
figure('units','centimeters','outerposition',[0 0 18 14]);
ax3=axes('Position',[.015 0.02 .5 1]);
hold on
thalweg_map(0,0,0,0);
m_text(-125.75,50.05,'a)','fontweight','bold','fontsize',8);
ax2=axes('Position',[.64 .58 .3 .38]);
hold on
ax1=axes('Position',[.64 .09 .3 .38]);
hold on
% 
% ax3=subplot(4,2,1:6);
% hold on
% thalweg_map(0,0,0,0);
% ax2=subplot(4,2,7);
% hold on
% ax1=subplot(4,2,8);
% hold on

axes(ax3);
for i=1:length(Dsort_min_day)
    m_scatter(Dsort_lon(i),Dsort_lat(i),25,mycol(i,:),'filled',...
        'markerfacealpha',0.75); 
end

axes(ax2)
for i=1:length(Dsort_min_day)
    plot(1:365,Dsort_harms(:,i),'linewidth',1,'color',mycol(i,:));
end

axes(ax1)
for i=1:length(Dsort_min_day)
    yyaxis left
    scatter(distance_N(i)/1000,Dsort_min_day(i),15,lnes(1,:),'filled');
    yyaxis right
    scatter(distance_N(i)/1000,Dsort_oamp(i),15,lnes(2,:),'filled');
end

axes(ax2)
ylabel('Dissolved Oxygen (\mumol kg^{-1})','fontsize',8,'fontweight','bold');
xlabel('Year day','fontsize',8,'fontweight','bold');
axis tight
grid on
text(0.9,0.9,'b)','units','normalized','fontweight','bold','fontsize',8);

axes(ax1)
yyaxis left
coeffs=polyfit(distance_N/1000,Dsort_min_day,1);
fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
    fittedy,'-','linewidth',2,'color',lnes(1,:));
ylabel('Minimum Year Day','fontsize',8,'fontweight','bold');
grid on
    
yyaxis right
coeffs=polyfit(distance_N/1000,Dsort_oamp,1);
fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
    fittedy,'-','linewidth',2,'color',lnes(2,:));
ylabel('Amplitude (\mumol kg^{-1})','fontsize',8,'fontweight','bold');    
xlabel('Northward Excursion (km)','fontsize',8,'fontweight','bold');
axis tight
text(0.1,0.9,'c)','units','normalized','fontweight','bold','fontsize',8);

%%
set(gcf, 'Color', 'w');
export_fig /ocean/sstevens/IW_project/figures/paper/ox_fits.png -png -m3

% In general there is a trend of decreasing amplitude and increasing
% yearday with northward excursion up the strait, however the relationship
% is much stronger in amplitude than yearday. This suggests that there is
% more across strait variability in yearday that amplitude, which is also
% seen in the kriging maps.

%% Kriging of oxmean
lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -122.500002324];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(station_lon,station_lat);

% dday = variogram([X' Y'],oxmean','nrbins',50,...
%     'maxdist',5e3,'plotit',false);
dday = variogram([X' Y'],oxmean');

[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,station_lon',station_lat',oxmean',grid_lon,grid_lat);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25;x],...
    [station_lat';49.833;y],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

celld_label=strrep(celld,'_ctd','');
celld_label=strrep(celld_label,'Noos','Nanoose');

figure('units','centimeters','outerposition',[0 0 11.4 12]);
% ax1=axes('Position',[0.3 0 0.4 1]);
ax1=axes('Position',[0 0 1 1]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(grid_lon,grid_lat,zi);
colormap(ax1,flipud(m_colmap('jet')));
m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');

[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,RGB);
end

markers='^dosh';
for i=1:length(celld)
    [X,Y]=m_ll2xy(station_lon(strcmp(rel_dataset,celld{i})),...
        station_lat(strcmp(rel_dataset,celld{i})));
    s(i)=scatter(X,Y,20,rgb('black'),markers(i),'filled');
end

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[ax,h]=m_contfbar(ax1,-0.03,[0.1 0.8],CS,CH,'endpiece','no','axfrac',.01,'fontsize',8);
ylabel(ax, {'DO mean/';'umol kg^{-1}'},'fontsize',8,'color','k','fontweight','bold');

bax=axes('Position',[0.367 0.15 .14 0.26]);
xticks([]);yticks([]);
box on
ax1=axes('Position',[0.4 0.2 .1 0.2]);
hold on
coeffs=polyfit(distance_N/1000,Dsort_oxmean,1);
scatter(distance_N/1000,Dsort_oxmean,30,'w','filled','markeredgecolor','r');
fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
    fittedy,'--','linewidth',3.5,'color',lnes(1,:));
ylabel('DO mean/ umol kg^{-1}','fontsize',8,'fontweight','bold');
xlabel('Distance North','fontsize',8,'fontweight','bold');
axis tight

%% Damping to represent mixing
% Create mixing 
mix=oamp./max(oamp);
Dsort_mix=mix(Nidx);

% Plot mix
dday = variogram([station_lon' station_lat'],mix','nrbins',50,...
    'maxdist',1,'plotit',false);

[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,0.4,400,dday.num,'plotit',false);

[zi,zivar] = kriging_use(vstruct,station_lon',station_lat',mix',grid_lon,grid_lat);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25],[station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

lat_lim=[48.5 50.3];
lon_lim=[-125.280000012121 -123.000002324];

celld_label=strrep(celld,'_ctd','');
celld_label=strrep(celld_label,'Noos','Nanoose');

figure('units','normalized','outerposition',[0 0 1 1]);
ax10=axes('Position',[0.3 0 0.4 1]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

Cm=flipud(cbrewer('div','RdYlBu',10));
[CS,CH]=m_contourf(grid_lon,grid_lat,zi);
colormap(Cm);
m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');

[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,RGB);
end

markers='^dosh';
for i=1:length(celld)
    [X,Y]=m_ll2xy(station_lon(strcmp(rel_dataset,celld{i})),...
        station_lat(strcmp(rel_dataset,celld{i})));
    s(i)=scatter(X,Y,20,rgb('black'),markers(i),'filled');
end

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[ax,h]=m_contfbar(ax10,-0.03,[0.1 0.8],CS,CH,'endpiece','no','axfrac',.01,'fontsize',8);
ylabel(ax, {'Oxygen Mixing'},'fontsize',8,'color','k','fontweight','bold');

bax=axes('Position',[0.367 0.15 .14 0.26]);
xticks([]);yticks([]);
box on
ax2=axes('Position',[0.4 0.2 .1 0.2]);
hold on
coeffs=polyfit(distance_N/1000,Dsort_mix,1);
scatter(distance_N/1000,Dsort_mix,30,'w','filled','markeredgecolor','r');
fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
    fittedy,'--','linewidth',3.5,'color',lnes(1,:));
ylabel('Oxygen Mixing','fontsize',8,'fontweight','bold');
xlabel('Distance North','fontsize',8,'fontweight','bold');
axis tight

%% Oxygen utilisation rate to gauge age
% Create delta O: the maximum mean oxygen (from the southeast end of the
% strait) minus the mean oxygen at each station.
deltaO=max(oxmean)-oxmean;

% Create age based on URate from Pawlowicz et al 2007, 0.9 ml/l/yr
% converted to umol/l/day ;
age=deltaO/((0.5*44.661)/365);
[~,idx]=max(min_day);
disp(['Water age from inflow to NW Tex.= ',num2str(age(idx))]);

Dsort_age=age(Nidx);

% % Plot AGE
% dday = variogram([station_lon' station_lat'],age','nrbins',50,...
%     'maxdist',1,'plotit',false);
% 
% [a,c,n,vstruct] = variogramfit(dday.distance,dday.val,0.4,400,dday.num,'plotit',false);
% 
% [zi,zivar] = kriging_use(vstruct,station_lon',station_lat',age',grid_lon,grid_lat);
% 
% % Define boundary of stations and NaN everything outside
% station_bound=alphaShape([station_lon';-125.25],[station_lat';49.833],1,'HoleThreshold',15);
% inside_idx=inShape(station_bound,grid_lon,grid_lat);
% zi(~inside_idx)=NaN;
% zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK
% 
% lat_lim=[48.5 50.3];
% lon_lim=[-125.280000012121 -123.000002324];
% 
% celld_label=strrep(celld,'_ctd','');
% celld_label=strrep(celld_label,'Noos','Nanoose');
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% ax1=axes('Position',[0.3 0 0.4 1]);
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
% hold on
% RGB=rgb('light grey');
% 
% Cm=flipud(cbrewer('div','RdYlBu',10));
% [CS,CH]=m_contourf(grid_lon,grid_lat,zi,[0:50:500]);
% colormap(Cm);
% m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');
% 
% [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
% k=[find(isnan(X(:,1)))];
% for i=1:length(k)-1
%     x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
%     y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
%     m_patch(x,y,RGB);
% end
% 
% markers='^dosh';
% for i=1:length(celld)
%     [X,Y]=m_ll2xy(station_lon(strcmp(rel_dataset,celld{i})),...
%         station_lat(strcmp(rel_dataset,celld{i})));
%     s(i)=scatter(X,Y,20,rgb('black'),markers(i),'filled');
% end
% 
% m_grid('linestyle','none','linewidth',2,'tickdir','out',...
%     'xaxisloc','bottom','yaxisloc','right','fontsize',14,'box','fancy');
% [ax,h]=m_contfbar(ax1,-0.03,[0.1 0.8],CS,CH,'endpiece','no','axfrac',.01,'fontsize',14);
% ylabel(ax, {'Age/ days'},'fontsize',14,'color','k','fontweight','bold');
% 
% bax=axes('Position',[0.367 0.15 .14 0.26]);
% xticks([]);yticks([]);
% box on
% ax1=axes('Position',[0.4 0.2 .1 0.2]);
% hold on
% coeffs=polyfit(distance_N/1000,Dsort_age,1);
% scatter(distance_N/1000,Dsort_age,30,'w','filled','markeredgecolor','r');
% fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
% l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
%     fittedy,'--','linewidth',3.5,'color',lnes(1,:));
% ylabel('Age/ days','fontsize',12,'fontweight','bold');
% xlabel('Distance North','fontsize',12,'fontweight','bold');
% axis tight
% 
% clc
% fprintf('Average along-strait IW speed= %2.1f cm/s,\n',coeffs(1)/(60*60*24)*100000);
% 
% % Assumes the utilisation rate is the same everywhere

