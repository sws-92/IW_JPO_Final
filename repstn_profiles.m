%% Plot IW details (seasonality, rep stations)
clc
clear
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
seas=load('/ocean/sstevens/IW_project/data/thick.mat');
seas.CD_profile(:,[15 132 139])=NaN; % Remove bad profile
station=load('sinfit_temp_data.mat');
Tstn=load('/ocean/sstevens/IW_project/data/represent_stns_t.mat');
idx=contains(station.rel_station_names,{'GO-4';'PR-2'});
station.station_lon(idx)=[];station.station_lat(idx)=[];
station.rel_dataset(idx)=[];
load BCcoast.mat

%% Plot seasonality
figure('units','centimeters','outerposition',[0 0 18 10.25],'color','w');

% Plot repstn
ax2=axes('Position',[0.075 0.11 .3 0.8]);
hold on
plot(1:365,Tstn.GO4.fitted_harm,'color','b','linewidth',1.5);
scatter(Tstn.GO4.day,Tstn.GO4.temp,30,'b','filled',...
    'MarkerEdgeColor','k','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.15);
plot(1:365,Tstn.PR2.fitted_harm,'color',rgb('dark green'),'linewidth',1.5);
scatter(Tstn.PR2.day,Tstn.PR2.temp,30,'filled','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.15,'markerfacecolor',rgb('dark green'));
errorbar(370,mean(Tstn.GO4.fitted_harm),range(Tstn.GO4.fitted_harm)/2,...
    'color','b')
scatter(370,mean(Tstn.GO4.fitted_harm),40,'b','marker','x');
errorbar(380,mean(Tstn.PR2.fitted_harm),range(Tstn.PR2.fitted_harm)/2,...
    'color',rgb('dark green'))
scatter(380,mean(Tstn.PR2.fitted_harm),40,'markeredgecolor',rgb('dark green'),'marker','x');
line([Tstn.GO4.coldday Tstn.GO4.coldday],[7.5 min(Tstn.GO4.fitted_harm)],...
    'color','b','linestyle','--','linewidth',2);
line([Tstn.PR2.coldday Tstn.PR2.coldday],[7.5 min(Tstn.PR2.fitted_harm)],...
    'color',rgb('dark green'),'linestyle','--','linewidth',2);
ylabel('Temperature (^oC)','fontsize',8,'fontweight','bold');
xlabel('Yearday','fontsize',8,'fontweight','bold');
axis tight
grid on
box on
text(0.9,0.05,'a)','units','normalized','fontweight','bold','fontsize',8);
xlim([0 400]);

% Plot quick map
ax2=axes('Position',[0.1 0.6 .075 0.33]);
lat_lim=[50.1 48.5];
    lon_lim=[-125.5 -122.522345235];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.25) 
hold on
% [X,Y]=m_ll2xy(ncst(:,1),ncst(:,2),'clip','patch');
% k=[find(isnan(X(:,1)))];
% for i=1:length(k)-1
%     x=ncst([k(i)+1:(k(i+1)-1) k(i)+1],1);
%     y=ncst([k(i)+1:(k(i+1)-1) k(i)+1],2);
%     m_patch(x,y,rgb('very light grey'));
% end
m_gshhs_f('patch',rgb('light grey'),'edgecolor','none');
m_northarrow(-125.32,49.75,0.15,'type',2,'linewidth',.25);
set(findall(gcf,'-property','FontSize'),'FontSize',8);

m_grid('linestyle','none','linewidth',0.25,'linecolor',rgb('light grey'),...
    'tickdir','out','xaxisloc','bottom','yaxisloc','left','fontsize',6,...
    'xtick',-126:-122,'ytick',48:51);
% Plot representative stations
a(2)=m_scatter(Tstn.GO4.lon,Tstn.GO4.lat,30,'b',...
    'filled','MarkerEdgeColor','k','markerfacealpha',0.7);
a(1)=m_scatter(Tstn.PR2.lon,Tstn.PR2.lat,30,'filled','MarkerEdgeColor',...
    'k','markerfacealpha',0.7,'markerfacecolor',rgb('dark green'));
% text(0.75,0.05,'b)','units','normalized','fontweight','bold','fontsize',8);


% Plot profiles
% coldday
dep=[20:250]';
% seas.CD_profile([1:19 201:400],:)=[];
% for i=1:size(seas.CD_profile,2)
%     idx=~isnan(seas.CD_profile(:,i));
%     sdm(:,i)=smooth(seas.CD_profile(idx,i),5);
% end
sdm=seas.CD_profile(20:250,:);
CD_dm=((sdm-mean(sdm,'omitnan'))./std(sdm,'omitnan'));
sm=smooth(mean(CD_dm,2,'omitnan'),5);
cm=linspecer(3);

SGidx=strcmp(seas.Dsort_dataset,'Stratogem_ctd');

ax1=axes('Position',[0.45 0.11 .25 0.8],'color','none');
plot(CD_dm,dep,'color',rgb('very light grey'));
hold on
plot(sm,dep,'w','linewidth',3)
plot(sm(1:39),dep(1:39),'k',sm(61:69),dep(61:69),'k',sm(61:69),dep(61:69),'k',...
    sm(91:99),dep(91:99),'k',sm(121:end),dep(121:end),'k','linewidth',1);
l(1)=plot(sm(40:60),dep(40:60),'color',cm(1,:),'linewidth',2);
l(2)=plot(sm(70:90),dep(70:90),'color',cm(2,:),'linewidth',2);
l(3)=plot(sm(100:120),dep(100:120),'color',rgb_x('goldenrod'),'linewidth',2);
% plot(smooth(mean(CD_dm(:,SGidx),2,'omitnan'),5),dep,'k:','linewidth',1)

set(ax1,'ydir','reverse','fontsize',8,'Xgrid','on'); 
axis tight;
ylim([20 250]);
xlim([-3 4]);
xlabel('Normalized Coldest Day','fontsize',8,'fontweight','bold');
ylabel('Depth (m)','fontsize',8,'fontweight','bold');
% xlim([-0.5 1]);
text(0.9,0.05,'b)','units','normalized','fontweight','bold','fontsize',8);
box off
set(ax1,'Color','none');

% amp
dep=[20:250]';
sdm=seas.tamp_profile(20:250,:);
tamp_dm=((sdm-mean(sdm,'omitnan'))./std(sdm,'omitnan'));
sm=smooth(mean(tamp_dm,2,'omitnan'),5);
cm=linspecer(3);

ax2=axes('Position',[0.725 0.11 .25 0.8]);
plot(tamp_dm,dep,'color',rgb('very light grey'));
hold on
l(1)=plot(sm,dep,'w','linewidth',3);
plot(sm(1:39),dep(1:39),'k',sm(61:69),dep(61:69),'k',sm(61:69),dep(61:69),'k',...
    sm(91:99),dep(91:99),'k',sm(121:end),dep(121:end),'k','linewidth',1);
l(1)=plot(sm(40:60),dep(40:60),'color',cm(1,:),'linewidth',2);
l(2)=plot(sm(70:90),dep(70:90),'color',cm(2,:),'linewidth',2);
l(3)=plot(sm(100:120),dep(100:120),'color',rgb_x('goldenrod'),'linewidth',2);
% plot(smooth(mean(tamp_dm(:,SGidx),2,'omitnan'),10),dep,'k:','linewidth',1)
set(ax2,'ydir','reverse','fontsize',8,'Xgrid','on'); axis tight;
ylim([20 250]);
ylabel('Depth/ m','fontsize',8,'fontweight','bold');
xlabel('Normalized Amplitude (\circC)','fontsize',8,'fontweight','bold');
xlim([-2 4]);
text(0.9,0.05,'c)','units','normalized','fontweight','bold','fontsize',8);
ax2.YAxis.Visible='off';
box off
set(ax2,'Color','none');

iax=axes('Position',[0.45 0.11 .525 0.8]);
ylim([20 250]);xlim([0 1]);set(iax,'ydir','reverse');
 idx=dep>=60 & dep<=80 | dep>=90 & dep<=110 | dep>=120 & dep<=140;
p(1)=patch([0 1 1 0],[60 60 80 80],cm(1,:),'edgecolor','none','facealpha',0.1);
p(2)=patch([0 1 1 0],[90 90 110 110],cm(2,:),'edgecolor','none','facealpha',0.1);
p(3)=patch([0 1 1 0],[120 120 140 140],rgb_x('goldenrod'),'edgecolor','none','facealpha',0.1);
linkaxes([iax ax1],'y'); 
iax.YAxis.Visible='off';iax.XAxis.Visible='off';
uistack(iax,'bottom');
iax.Color='none'; iax.YGrid='on';

%%
export_fig /ocean/sstevens/IW_project/figures/paper/repstn_profiles.png -png -m3
export_fig /ocean/sstevens/IW_project/figures/paper/vec/repstn_profiles.pdf -dpdf
%% Contour plot 
didx=dep>40;
dep=dep(didx);
Dsort_ampdm=tamp_dm(didx,seas.Nidx);
Dsort_CDdm=CD_dm(didx,seas.Nidx);
X=0:15:max(seas.distance_N)/1000;

SGidx=strcmp(seas.Dsort_dataset,'Stratogem_ctd');
SGdis=mean(seas.distance_N(SGidx))/1000;

amp_cont=NaN(length(dep),length(X));

for i=1:length(X)-1
    idx=seas.distance_N/1000>X(i) & seas.distance_N/1000<X(i+1);
    amp_cont(:,i)=mean(Dsort_ampdm(:,idx),2,'omitnan');
end

Xind=find(X>SGdis);

figure('units','centimeters','outerposition',[0 0 13.5 11.5],'color','w');
% axes('Position',[0.15 0.57 0.10 0.09]);
% 
% axes('Position',[0.15 0.57 0.10 0.09]);
hold on
contourf(X,dep,amp_cont,[-2 -1:0.25:2],'LineColor','none');
caxis([-1 1]);
colormap(cmocean('balance','pivot',-0.5));
% colormap(cmocean('amp'));
set(gca,'ydir','reverse','fontsize',8);axis tight;
% [C,h]=contour(X,dep,amp_cont,[-0.4 -0.4],'color',...
%     [0.4 0.4 0.4],'linewidth',1);
xlabel('Northward Excursion (km)','Fontweight','bold');
ylabel('Depth (m)','fontweight','bold');
xl=ylim;
plot([X(Xind(1)) X(Xind(1))],[min(xl) max(xl)],'color',[0.4 0.4 0.4],...
    'linestyle','--');


%% Contour plot 
Dsort_CDdm=CD_dm(didx,seas.Nidx);
X=0:20:max(seas.distance_N)/1000;

CD_cont=NaN(length(dep),length(X));

for i=1:length(X)-1
    idx=seas.distance_N/1000>X(i) & seas.distance_N/1000<X(i+1);
    CD_cont(:,i)=mean(Dsort_CDdm(:,idx),2,'omitnan');
end

figure('units','centimeters','outerposition',[0 0 13.5 11.5],'color','w');
% axes('Position',[0.15 0.57 0.10 0.09]);
% 
% axes('Position',[0.15 0.57 0.10 0.09]);
hold on
contourf(X,dep,CD_cont,[-2 -0.5:0.25:2],'LineColor','none');
caxis([-1 1]);
% colormap(cmocean('balance','pivot',-0.5));
colormap(cmocean('amp'));
set(gca,'ydir','reverse','fontsize',8);grid on; axis tight;
% [C,h]=contour(X,dep,amp_cont,[-0.4 -0.4],'color',...
%     [0.4 0.4 0.4],'linewidth',1);
xlabel('Northward Excursion (km)','Fontweight','bold');
ylabel('Depth (m)','fontweight','bold');
xl=ylim;
plot([SGdis SGdis],[min(xl) max(xl)],'color',[0.4 0.4 0.4],...
    'linestyle','--');