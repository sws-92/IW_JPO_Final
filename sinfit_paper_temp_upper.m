%% Fit sinusoids to all stations and use objective mapping to to create map
% of coldest yearday for entire strait

% LOOK AT TIMESCALE OF WATER TRAVELLING THROUGH BC COMPARED TO OUTSIDE BC

% Stations to use
% Stratogem
% PSF
% Hakai
% IOS
% any Nanoose?
clc
clear
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
load /ocean/sstevens/IW_project/data/all_ctd_July.mat
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
%     all_ctd.temp(:,idx)=[];
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
warning('off','last')

%% Create sinusoid fits to find coldest yearday

% Loop through all stations with more than 10 profiles
cold_day=NaN(1,length(rel_station_names));
fitted_harms=NaN(365,length(rel_station_names));
station_lat=NaN(1,length(length(rel_station_names)));
station_lon=NaN(1,length(length(rel_station_names)));
bad_switch=NaN(1,length(length(rel_station_names)));
h=NaN(1,length(length(rel_station_names)));
pval=NaN(1,length(length(rel_station_names)));
amp=NaN(1,length(rel_station_names));
cold_dayCI=NaN(length(rel_station_names),2);
ampCI=NaN(length(rel_station_names),2);
cold_daySD=NaN(length(rel_station_names),1);
ampSD=NaN(length(rel_station_names),1);
HS_date=[];
HS_depth=[];
HS_coords=[];  


% % q=input('Do you want to plot stations and sinusoids? y/n: ','s');
% % if strcmp(q,'y')
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     ax3=subplot(3,4,[6:8 10:12]);
%     thalweg_map(0,0,0,0);
%     ax2=subplot(3,4,1:4);
%     ax1=subplot(3,4,[5 9]);
%     set(gca,'ydir','reverse');
% % end
% q='n';

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
    nmsk=logical(sum(isnan(all_ctd.temp(90:110,stn_idx)))>3);
    all_ctd.temp(:,nmsk)=NaN;
    
    tmp_temp=nanmean(all_ctd.temp(60:80,stn_idx));
    tmp_temp(nan_msk)=[];
    tmp_temp=tmp_temp(sidx);
    tmp_temp=inpaint_nans(tmp_temp');
    
    % fit sinusoid
    [amp(i),phase,frac,offset,yy]=fit_harmonics(tmp_temp,tmp_day',...
        1,365,0.01);
    
    for tt=1:365
        fitted_harm(tt)=offset+amp(i)*cos(2*pi*tt/365 + phase);
    end
    
    [~,cold_day(i)]=min(fitted_harm);
    fitted_harms(:,i)=fitted_harm;
    
    % Check goodness-of-fit
    h(i)=nanmean(abs(tmp_temp-interp1(1:365,fitted_harm,tmp_day)));
    
    % Bootstrapping
    % Create sample indices
    I=bootrnd(length(tmp_temp),500);
    BSamp=NaN(length(tmp_temp),1);
    BSphase=NaN(length(tmp_temp),1);
    BScold_day=NaN(length(tmp_temp),1);
    
    % Loop through bootstrap samples and do harmonic analysis
    for ii=1:size(I,2)
        sample=sortrows([tmp_temp(I(:,ii)) tmp_day(I(:,ii))],2);
        [BSamp(ii),BSphase]=fit_harmonics(sample(:,1),sample(:,2),...
            1,365,0.01);
        for tt=1:365
            BSfitted_harm(tt)=offset+BSamp(ii)*cos(2*pi*tt/365 + BSphase);
        end
        [~,BScold_day(ii)]=min(BSfitted_harm);
    end
    
    % Find bootstrap 95% CI
    SEM=std(BSamp)/sqrt(length(BSamp));
    ts = tinv([0.025  0.975],length(BSamp)-1);
    ampCI(i,:)=nanmean(BSamp)+ts*SEM;
    ampSD(i)=nanstd(BSamp);
    
    SEM=std(BScold_day)/sqrt(length(BScold_day));
    ts = tinv([0.025  0.975],length(BScold_day)-1);
    cold_dayCI(i,:)=nanmean(BScold_day)+ts*SEM;
    cold_daySD(i)=nanstd(BScold_day);
    
%     if station_lat(i)<49 %strcmp(q,'y') && cold_day(i)~=1
% %         axes(ax3);
% %         s=m_scatter(station_lon(i),station_lat(i),50,'filled');
% %         axes(ax2)
% %         p(1)=plot(1:365,fitted_harm);
% %         hold on
% %         p(2)=scatter(tmp_day,tmp_temp);
% %         p(3)=scatter(cold_day(i),fitted_harm(cold_day(i)),[],'b','filled');
% %         l=legend(rel_dataset{i});
% %         axes(ax1)
% %         plot(nanmean(all_ctd.temp(:,stn_idx),2),1:400);
% %         set(gca,'ydir','reverse');
% %         keyboard
% %         %         q2=input('Is this a good station? y/n: ','s');
% %         %         if strcmp(q2,'n')
% %         %             bad_switch(idx)=1;
% %         %         end
% %         delete(s)
% %         delete(p)
% %         delete(l)
%         HS_date=[HS_date;tmp_mtime'];
%         HS_depth=[HS_depth;...
%             find(isfinite((nanmean(all_ctd.temp(:,stn_idx),2)))==1,1,'last')];
%         HS_coords=[HS_coords; [station_lon(i) station_lat(i)]];
%     end
    
    % Save certain stations for plotting elsewhere
    if strcmp(rel_station_names(i),'PR-2')
        PR2.amp=amp(i); PR2.fitted_harm=fitted_harm;
        PR2.day=tmp_day; PR2.temp=tmp_temp;
        PR2.lat=station_lat(i);PR2.lon=station_lon(i);
        PR2.coldday=cold_day(i); keyboard
        
    elseif strcmp(rel_station_names(i),'GO-4')
        GO4.amp=amp(i); GO4.fitted_harm=fitted_harm;
        GO4.day=tmp_day; GO4.temp=tmp_temp;
        GO4.lat=station_lat(i);GO4.lon=station_lon(i);
        GO4.coldday=cold_day(i); keyboard
        
    elseif strcmp(rel_station_names(i),'S22')
        S22.amp=amp(i); S22.fitted_harm=fitted_harm;
        S22.day=tmp_day; S22.temp=tmp_temp;
        S22.lat=station_lat(i);S22.lon=station_lon(i);
        S22.coldday=cold_day(i);
    end
    

end

% % Remove broken stations
idx=cold_day==1 | amp>2;

% Remove Texada Island and Saanich stations
idx_TI_Saan=contains(rel_station_names,{'IOS_22';'LS-3';'LS-4';'LS-5';'BS-6'});
idx=logical(idx+idx_TI_Saan);

cold_day(idx)=[]; station_lat(idx)=[]; station_lon(idx)=[];
rel_dataset(idx)=[]; rel_station_names(idx)=[]; fitted_harms(:,idx)=[];
rel_stations(idx)=[]; amp(idx)=[]; ampSD(idx)=[]; cold_daySD(idx)=[]; 
ampCI(idx,:)=[]; cold_dayCI(idx,:)=[]; 

% Save data for plotting elsewhere
save('/ocean/sstevens/IW_project/data/represent_stns_t.mat','PR2','GO4','S22');
save('/ocean/sstevens/IW_project/data/sinfit_temp_data.mat');


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

distance_idx=NaN(1,length(station_lat));
for i=1:length(station_lat)
    [~,idxb]=min(abs(center_line(2,:)-station_lat(i)));
    distance_N(i)=center_line_dist(idxb);
end

% Sort by distance
[distance_N,Nidx]=sort(distance_N);
Dsort_coldday=cold_day(Nidx);
Dsort_amp=amp(Nidx);
Dsort_harms=fitted_harms(:,Nidx);
Dsort_lat=station_lat(Nidx);
Dsort_lon=station_lon(Nidx);
Dsort_ampSD=ampSD(Nidx);
Dsort_cold_daySD=cold_daySD(Nidx);
Dsort_dataset=rel_dataset(Nidx);
Dsort_stations=rel_station_names(Nidx);

save('/ocean/sstevens/IW_project/data/sinfit_temp_data_plotting.mat');

%% Plotting
pth=load('pathXY.mat');

% Set up XY projection
lat_limV=[min(grid_lat(:)) max(grid_lat(:))];
lon_limV=[min(grid_lon(:)) max(grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(station_lon,station_lat);
[Xgrid,Ygrid]=m_ll2xy(grid_lon,grid_lat);

dday = variogram([X' Y'],cold_day','maxdist',49000,'plotit',false);

[a,c,~,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',cold_day',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25],...
    [station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK
zi_age=zi-40;

lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -122.500002324];
 
celld_label=strrep(celld,'_ctd','');
celld_label=strrep(celld_label,'Noos','Nanoose');

figure('units','centimeters','outerposition',[0 0 16 22]);
ax3=axes('Position',[0.05 0.55 1 0.45]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(grid_lon,grid_lat,zi_age,[0:10:20 30:15:130],'linecolor','k');
colormap(m_colmap('jet'));
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
m_scatter(station_lon([114 32]),station_lat([114 32]),7,...
    'markerfacecolor','non','markeredgecolor','w');
m_text(-122.75,50,'a)','fontweight','bold','fontsize',8);

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[ax,h]=m_contfbar(ax3,0.11,[0.0 0.95],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
ylabel(ax,{'Age (days)'},'fontsize',8,'color','k','fontweight','bold');
m_plot(pth.Xlon,pth.Xlat,'w:','linewidth',1.5);

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

damp = variogram([X' Y'],amp',...
    'maxdist',40000,'plotit',false);

[a,c,n,vstruct] = variogramfit(damp.distance,damp.val,[],[],[],'plotit',false);

[zi,zivar] = kriging_use(vstruct,X',Y',amp',Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25;],...
    [station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK
zi_amp=1.5-zi;

% Set up figure
lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -122.500002324];

% Set up figure
ax1=axes('Position',[0.05 0.05 1 0.45]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(grid_lon,grid_lat,zi,[0:0.2:0.8 1:.25:1.5]); % FIX THIS
colormap(ax1,flipud(m_colmap('jet')));
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
m_text(-122.75,50,'c)','fontweight','bold','fontsize',8);

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[ax,h]=m_contfbar(ax1,0.11,[0.0 0.95],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
ylabel(ax,{'Amplitude/ ^oC'},'fontsize',8,'color','k','fontweight','bold');

% Plot semivariogram
axes('Position',[0.255 0.06 0.28 0.2]);
set(gca, 'box','on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[]);
text(0.1,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);
axes('Position',[0.32 0.095 0.2 0.15]);
variogramfit_plot(damp.distance/1000,damp.val,[],[],[],'plotit',true);
set(gca,'fontsize',8);

% caxis([0 1.5]);
% 
% axes(ax1)
% l=legend(s,celld_label,'location','southwest','fontsize',7);
% title(l,'Dataset','fontsize',8);


% for i=1:length(X)
%     i
%     res(i)=cold_day(i)-griddata(grid_lon,grid_lat,zi,station_lon(i),station_lat(i));
%     res(i)
% end

% load krig_res.mat
% msk=abs(res)>6;
% m_scatter(X(msk),Y(msk),100000,'w','filled');

% EBC advection scale
S_idx=strcmp(rel_station_names,'SN-1');
N_idx=strcmp(rel_station_names,'IOS_44');
EBCtime=cold_day(S_idx)-cold_day(N_idx); %days
EBCdist=sum(gsw_distance(pth.Xlon,pth.Xlat)/1000);%KM
EBCspeed=(EBCdist*100000)/(EBCtime*60*60*24); %cm s-1

fprintf('EBC distance= %2.0f km \n EBC time= %2.0f days \n EBC speed= %2.2f cm s-1 \n',...
    EBCdist,EBCtime,EBCspeed);

%%
set(gcf, 'Color', 'w');
export_fig /ocean/sstevens/IW_project/figures/paper/temp_colday_amp.png -png -m3

%% Plot sinusoids
clear global
% Find station distances along center line
distance_idx=NaN(1,length(station_lat));
for i=1:length(station_lat)
    [~,idxb]=min(abs(center_line(2,:)-station_lat(i)));
    distance_N(i)=center_line_dist(idxb);
end

% Sort by distance
[distance_N,Nidx]=sort(distance_N);
Dsort_coldday=cold_day(Nidx);
Dsort_amp=amp(Nidx);
Dsort_harms=fitted_harms(:,Nidx);
Dsort_lat=station_lat(Nidx);
Dsort_lon=station_lon(Nidx);
Dsort_ampSD=ampSD(Nidx);
Dsort_cold_daySD=cold_daySD(Nidx);
Dsort_dataset=rel_dataset(Nidx);
Dsort_stations=rel_station_names(Nidx);

% Longest timescale? Southernmost station to maximum (which is NW of
% Texada)
templong=max(cold_day)-Dsort_coldday(1);
disp(['Inflow to NW time (coldday)=',num2str(templong)]);

% create colmap
mycol=m_colmap('jet',length(Dsort_coldday));

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
for i=1:length(Dsort_coldday)
    m_scatter(Dsort_lon(i),Dsort_lat(i),25,mycol(i,:),'filled',...
        'markerfacealpha',0.75); 
end

axes(ax2)
for i=1:length(Dsort_coldday)
    plot(1:365,Dsort_harms(:,i),'linewidth',1,'color',mycol(i,:));
end

axes(ax1)
for i=1:length(Dsort_coldday)
    yyaxis left
    scatter(distance_N(i)/1000,Dsort_coldday(i),15,lnes(1,:),'filled');
    yyaxis right
    scatter(distance_N(i)/1000,Dsort_amp(i),15,lnes(2,:),'filled');
end

axes(ax2)
ylabel('Temperature (^oC)','fontsize',8,'fontweight','bold');
xlabel('Yearday','fontsize',8,'fontweight','bold');
axis tight
grid on
text(0.1,0.9,'b)','units','normalized','fontweight','bold','fontsize',8);

axes(ax1)
yyaxis left
coeffs=polyfit(distance_N/1000,Dsort_coldday,1);
fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
    fittedy,'-','linewidth',2,'color',lnes(1,:));
ylabel('Coldest Yearday','fontsize',8,'fontweight','bold');
grid on
    
yyaxis right
coeffs=polyfit(distance_N/1000,Dsort_amp,1);
fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
    fittedy,'-','linewidth',2,'color',lnes(2,:));
ylabel('Amplitude (^oC)','fontsize',8,'fontweight','bold');    
xlabel('Northward Excursion (km)','fontsize',8,'fontweight','bold');
axis tight
text(0.1,0.1,'c)','units','normalized','fontweight','bold','fontsize',8);


%%
set(gcf, 'Color', 'w');
% export_fig /ocean/sstevens/IW_project/figures/paper/temp_fits.eps -depsc
export_fig /ocean/sstevens/IW_project/figures/paper/temp_fits.png -png -m3
    
% In general there is a trend of decreasing amplitude and increasing
% yearday with northward excursion up the strait, however the relationship
% is much stronger in amplitude than yearday. This suggests that there is
% more across strait variability in yearday that amplitude, which is also
% seen in the kriging maps.

%% Plot bootstrapping error estimates

% Set up XY projection
lat_limV=[min(grid_lat(:)) max(grid_lat(:))];
lon_limV=[min(grid_lon(:)) max(grid_lon(:))];
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(station_lon,station_lat);
[Xgrid,Ygrid]=m_ll2xy(grid_lon,grid_lat);

dday = variogram([X' Y'],cold_daySD,'plotit',true);

[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',true);

[zi,zivar] = kriging_use(vstruct,X',Y',ampSD,Xgrid,Ygrid);

% Define boundary of stations and NaN everything outside
station_bound=alphaShape([station_lon';-125.25],[station_lat';49.833],1,'HoleThreshold',15);
inside_idx=inShape(station_bound,grid_lon,grid_lat);
zi(~inside_idx)=NaN;
zi(Zint>-10)=NaN; % 10m SHALLOW WATER MASK

lat_lim=[48.4 50.15];
lon_lim=[-125.280000012121 -123.000002324];

celld_label=strrep(celld,'_ctd','');
celld_label=strrep(celld_label,'Noos','Nanoose');

figure('units','normalized','outerposition',[0 0 1 1]);
ax3=axes('Position',[0.05 0 0.4 1]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(grid_lon,grid_lat,zi,0:0.05:0.3);
colormap(ax3,m_colmap('jet'));
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
[ax,h]=m_contfbar(ax3,-0.03,[0.1 0.8],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
ylabel(ax, {'Coldest';'yearday'},'fontsize',8,'color','k','fontweight','bold');
clc


%% Test confidence Kriging
m_proj('UTM','lon',lon_limV,'lat',lat_limV);   % Projection
[X,Y]=m_ll2xy(station_lon,station_lat);
[Xgrid,Ygrid]=m_ll2xy(grid_lon,grid_lat);
dday = variogram([X' Y'],cold_day','nrbins',50,...
    'maxdist',1e5,'plotit',false);
[a,c,n,vstruct] = variogramfit(dday.distance,dday.val,[],[],[],'plotit',false);
[zi,zivar] = kriging_use(vstruct,X',Y',cold_day',Xgrid,Ygrid);
lon_vec=Xgrid(:);
lat_vec=Ygrid(:);

for i=1:length(X)
    i
%   waitbar(h,i/length(X),[num2str(i),'/',num2str(length(X))]);
    [~,idx1]=min(abs(lon_vec-X(i)));
    [~,idx2]=min(abs(lat_vec-Y(i)));
    
    res(i)=cold_day(i)-griddata(Xgrid,Ygrid,zi,lon_vec(idx1),lat_vec(idx2));
    res(i)
end

load krig_res.mat

figure('units','normalized','outerposition',[0 0 1 1]);
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');
% CH=m_pcolor(grid_lon,grid_lat,zi,'levels',40:15:140,'facealpha',0.5);
% m_contour(grid_lon,grid_lat,Zint,[-10 -10],'k');
[C,h]=m_contour(grid_lon,grid_lat,zi,40:15:140);
colormap(m_colmap('jet'));
clabel(C,h,'fontsize',8);
for i=1:6
    sbehind(i)=plot(-1,-1,'o','markersize',i*15*exp(-i/6)*0.8,...
        'color','k','markerfacecolor','k');
    err{i}=['Error=',num2str(i),' days'];
end

[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,RGB);
end

for i=1:length(res)
    [X,Y]=m_ll2xy(station_lon(i),station_lat(i));
    s(i)=plot(X,Y,'o','markersize',abs(res(i))*15*exp(-abs(res(i))/6)*0.8,'color','k',...
        'markerfacecolor','k');
end

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[res_sort,sidx]=sort(abs(res));
s=s(sidx);

l=legend(sbehind,err,'location','southwest');
title(l,'Coldest day');
siz=l.Position;
set(l,'position',[siz(1)*1.05 siz(2)*1.1 siz(3) siz(4)*2.25]);

%% Damping to represent mixing
% Create mixing 
mix=amp./max(amp);
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

% 
% lat_lim=[48.5 50.15];
% lon_lim=[-125.280000012121 -123.000002324];
% 
% celld_label=strrep(celld,'_ctd','');
% celld_label=strrep(celld_label,'Noos','Nanoose');
% 
% figure('units','normalized','outerposition',[0 0 1 1]);
% ax3=axes('Position',[0.3 0 0.4 1]);
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
% hold on
% RGB=rgb('light grey');
% 
% Cm=flipud(cbrewer('div','RdYlBu',10));
% [CS,CH]=m_contourf(grid_lon,grid_lat,zi);
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
%     'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
% [ax,h]=m_contfbar(-0.03,[0.1 0.8],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
% ylabel(ax, {'Temperature Mixing'},'fontsize',8,'color','k','fontweight','bold');
% 
% bax=axes('Position',[0.367 0.15 .14 0.26]);
% xticks([]);yticks([]);
% box on
% ax1=axes('Position',[0.4 0.2 .1 0.2]);
% hold on
% coeffs=polyfit(distance_N/1000,Dsort_mix,1);
% scatter(distance_N/1000,Dsort_mix,30,'w','filled','markeredgecolor','r');
% fittedy=polyval(coeffs,linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)));
% l1=plot(linspace(min(distance_N/1000),max(distance_N/1000),length(distance_N)),...
%     fittedy,'--','linewidth',3.5,'color',lnes(1,:));
% ylabel('Temperature Mixing','fontsize',8,'fontweight','bold');
% xlabel('Distance North','fontsize',8,'fontweight','bold');
% axis tight
% [r,p]=corrcoef(fittedy,Dsort_mix')


%% Plot histograms of southern stations
HS_date_hist=datevec(sort(HS_date));

figure;
subplot(5,1,1);
hist(HS_date_hist(:,2),12);
title('Southern Stations (<49\circN) sampling');
xlabel('Month');
ylabel('Profiles');

subplot(5,1,2);
hist(HS_date_hist(:,1));
xlabel('Years');

subplot(5,1,3);
hist(HS_depth,20);
xlabel('Max depth/ m');

subplot(5,1,4:5);
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
load /ocean/sstevens/IW_project/data/thalweg.mat
fname='/ocean/rich/more/mmapbase/noaa_bc3/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');
lon=ncread(fname,'lon');
lat_lim=[48.3 49.1];
lon_lim=[-123.569298287 -122.22345235];
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);

Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);
load('BCcoast');

% Set up figure
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on

caxis([-400 0]);
[C,h]=m_contour(lon(ilon),lat(ilat),Z',0:50:300,'color',rgb('light grey'));
clabel(C,h);
[X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
k=[find(isnan(X(:,1)))];
for i=1:length(k)-1
    x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
    y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
    m_patch(x,y,rgb('light grey'));
end
m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');

m_scatter(HS_coords(:,1),HS_coords(:,2),50,'k','marker','x');


%% Plot rate of amplitude loss
amp_loss=zi_amp./zi_age;
amp_loss(amp_loss<0 | amp_loss>0.03)=NaN;

figure('units','centimeters','outerposition',[0 0 17 25]);
ax3=axes;
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
hold on
RGB=rgb('light grey');

[CS,CH]=m_contourf(grid_lon,grid_lat,amp_loss,0:0.003:0.03,'linecolor','k');
colormap(cmocean('amp'));
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

m_grid('linestyle','none','linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',8,'box','fancy');
[ax,h]=m_contfbar(ax3,0.11,[0.0 0.95],CS,CH,'endpiece','no','axfrac',.03,'fontsize',8);
ylabel(ax,{'Coldest yearday'},'fontsize',8,'color','k','fontweight','bold');



