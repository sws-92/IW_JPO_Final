
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));

%% Load bathy data
clear
fname='/ocean/rich/more/mmapbase/noaa_bc3/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');
lon=ncread(fname,'lon');
lat_lim=[48 50.3];
lon_lim=[-125.99298287 -122.22345235];
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z_lon=lon(ilon);
Z_lat=lat(ilat);
% [Z_lon, Z_lat]=meshgrid(Z_lon,Z_lat);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

% Find bathy slope direction
spac=Z_lon(2)-Z_lon(1);
[aspect,~,gradN,gradE]=gradientm(Z,[1/spac max(Z_lat(:)) min(Z_lon(:))]);

[Z_lon, Z_lat]=meshgrid(Z_lon,Z_lat);

% % Check aspect
% m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
% figure;
% m_gshhs_l('patch',rgb('light grey'),'edgecolor','k');
% hold on;
% m_contourf(Z_lon,Z_lat,fliplr(rot90(gradE,3)),'linestyle','none');
% m_grid('linestyle','none','linewidth',0.1,'tickdir','out');

aspect=fliplr(rot90(aspect,3));
gradN=fliplr(rot90(gradN,3));
gradE=fliplr(rot90(gradE,3));

% Convert to 0-90?

clearvars -except Z_lon Z_lat Z aspect gradN gradE
%%
% load DP_2016age.mat
dir_list=dir('/ocean/sstevens/ariane_old/results_BP/2016*365d');

lon_lim=[-127 -121];
lat_lim=[46.5 52];
clear global

m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection

% make axis and do rotation
yl=47.85:0.01:51;
xl=ones(1,length(yl))*-124.0342;
vl=[xl;yl];
yl_cent=yl(round(length(yl)/2));xl_cent=xl(round(length(yl)/2));
hl=[-125.4:0.01:-122.4 ; ones(1,length(-125.4:0.01:-122.4))*yl_cent];
lcenter=repmat([xl_cent; yl_cent],1,length(xl));
R=[cosd(60) -sind(60);sind(60) cosd(60)];
sl=vl-lcenter;lso=R*sl;AS_axis=lso+lcenter;
lcenter=repmat([xl_cent; yl_cent],1,length(hl));
sl=hl-lcenter;lso=R*sl;AcS_axis=lso+lcenter;

[x,y]=m_ll2xy(xl,yl);
v=[x;y];
y_cent=y(round(length(y)/2));x_cent=x(round(length(y)/2));
lnes=lines;

middle=[x_cent*6370e3; y_cent*6370e3];

% % Check axis
% figure;
% thalweg_map(0,0,0,0);
% m_plot(x,y);hold on,m_plot(AS_axis(1,:),AS_axis(2,:));
% m_scatter(x_cent,y_cent,20,'filled');

load PT_AS_dist.mat
gate=[X(1:2) Y(1:2)];
station_bound=alphaShape(X,Y,1,'HoleThreshold',15);

tic

% % Map test  UNCOMMENT FOR TEST
% figure;
% thalweg_map(0,0,0,0);
% m_plot(AS_axis(1,:),AS_axis(2,:),'linewidth',1,'color',[0 0 0],'linest','--');
% m_plot(gate(:,1),gate(:,2),'k','linewidth',2);

% Rotation axis for trajectories
R=[cosd(-60) -sind(-60);sind(-60) cosd(-60)];

res=1/24; % Hourly res
ydisp=single(NaN(length(res:res:365),20352));
xdisp=single(NaN(length(res:res:365),20352));
% across_disp=NaN(length(res:res:365),20352);
% along_disp=NaN(length(res:res:365),20352);

load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;

count2=0;
for i=1:length(dir_list)
    
    % Load data
    fpath = [dir_list(i).folder '/'  dir_list(i).name '/'];
    filestr=[fpath 'TrajTBL'];
    load(filestr);
    disp(i)
    toc
    
    for j=1:length(TrajTBL)
        count=0;
        mn_traj=NaN(length(res:res:365),3);
        
        for ii=res:res:365
            count=count+1;
            mask=(TrajTBL(j).Age_days>ii-res & TrajTBL(j).Age_days<ii);
            mn_traj(count,1)=mean(TrajTBL(j).lon(mask));
            mn_traj(count,2)=mean(TrajTBL(j).lat(mask));
            mn_traj(count,3)=mean(TrajTBL(j).z_m(mask))*-1;
        end
        qtraj=mn_traj;
        
        msk=isnan(mn_traj(:,1));
        if sum(msk)
            mn_traj(msk,:)=0;
        end
        
        msk=mn_traj(:,3)<50 | mn_traj(:,3)>150;
        if sum(msk)
            mn_traj(msk,:)=0;
        end
        
        msk=~inShape(station_bound,mn_traj(:,1),mn_traj(:,2));
        mn_traj(msk,:)=NaN;
        
%         figure;
%         m_gshhs_l('patch',rgb('light grey'),'edgecolor','k');
%         hold on;
%         m_plot(qtraj(:,1),qtraj(:,2),'color',[0.8 0.8 0.8]);
%         [x,y]=m_ll2xy(AS.lon,AS.lat);
%         station_bound=alphaShape(x,y);
%         station_bound.Alpha=station_bound.Alpha*2;
%         plot(station_bound,'FaceColor','red','FaceAlpha',0.1,'edgecolor','none');
%         m_plot(qtraj(~msk,1),qtraj(~msk,2),'r');
        
        XYtraj=[];
        [XYtraj(1,:),XYtraj(2,:)]=m_ll2xy(mn_traj(:,1),mn_traj(:,2));
        
        % Do rotation
        center=repmat([x_cent; y_cent],1,length(mn_traj));
        sa=XYtraj-center;so=R*sa;AS_traj=so+center;
       
        count2=count2+1;
        xdisp(:,count2)=single((AS_traj(1,:)-AS_traj(1,1))*6370e3);
        ydisp(:,count2)=single((AS_traj(2,:)-AS_traj(2,1))*6370e3);
        
%         % Find displacement relative to depth contours
%         a=diff(XYtraj(1,:)-XYtraj(1,1))*6370e3; % displacement E in m per t relative to t=0
%         b=diff(XYtraj(2,:)*6370e3-XYtraj(2,1)*6370e3); % displacement N per t relative to t=0
%         
%         for jj=2:length(mn_traj)
%             %             aspect_i=qinterp2(Z_lon,Z_lat,aspect,mn_traj(jj,1),...
%             %                 mn_traj(jj,2),0)*-1;
%             %
%             %             % Do rotation
%             %             Ro=[cosd(aspect_i) -sind(aspect_i);sind(aspect_i) cosd(aspect_i)];
%             %             sa=[a(jj-1);b(jj-1)]; so=Ro*sa;
%             
%             %             across_disp(jj,count2)=so(1);
%             %             along_disp(jj,count2)=so(2);
%             
%             c=qinterp2(Z_lon,Z_lat,gradE,mn_traj(jj,1),...
%                 mn_traj(jj,2),0);
%             d=qinterp2(Z_lon,Z_lat,gradN,mn_traj(jj,1),...
%                 mn_traj(jj,2),0);
%             
%             e=fliplr([a(jj-1) b(jj-1);c/norm(c) d/norm(d)]);
% %             e=[a(jj-1) b(jj-1);c/norm(c) d/norm(d)];
%             
%             along_disp(jj,count2)=det(e);
%             across_disp(jj,count2)=dot(e(1,:),e(2,:));
%         end
        
%         along_disp(:,count2)=cumsum(along_disp(:,count2),'omitnan');
%         across_disp(:,count2)=cumsum(across_disp(:,count2),'omitnan');
%         along_disp(msk,count2)=NaN;across_disp(msk,count2)=NaN;
%         
%         
% %        eyeball the rotation
%         figure;
%         m_gshhs_l('patch',rgb('light grey'),'edgecolor','k');
%         hold on;
%         m_plot(qtraj(:,1),qtraj(:,2),'color',[0.8 0.8 0.8]);
%         m_plot(AS_axis(1,:),AS_axis(2,:));
%         m_plot(vl(1,:),vl(2,:));
%         scatter(x_cent,y_cent,20,'filled');
%         plot(XYtraj(1,:),XYtraj(2,:));
%         plot(AS_traj(1,:),AS_traj(2,:));
%         m_grid('linestyle','none','linewidth',0.1,'tickdir','out',...
%             'xticklabels',[],'yticklabels',[]);
%         
%         figure
%         subplot(1,2,1)
%         plot(across_disp(:,count2),'color',lnes(1,:));
%         hold on
%         plot(xdisp(:,count2),'color',lnes(2,:));
%         plot(along_disp(:,count2),'color',lnes(3,:));
%         plot(ydisp(:,count2),'color',lnes(4,:));
%         grid on;
%         legend('across','x','along','y','location','best');
%         
%         subplot(1,2,2)
%         scatter(xdisp(:,count2),ydisp(:,count2),5,'filled');
%         hold on
%         scatter(xdisp(1,count2),ydisp(1,count2),10,'r','filled');
%         grid on;ylabel('y displacement (m)');xlabel('x displacement (m)');
%         keyboard
%         close all

    end
end

toc
clearvars TrajTBL 

xvel=diff(xdisp); % m/hr
yvel=diff(ydisp); % m/hr
% across_vel=diff(across_disp); % m/hr
% along_vel=diff(along_disp);

save('/ocean/sstevens/ariane_old/results_BP/gate_cross_v3.mat');

%% Plot stats
clear
load /ocean/sstevens/ariane_old/results_BP/gate_cross_v3.mat
xdisp=double(xdisp); xvel=double(xvel);
ydisp=double(ydisp); yvel=double(yvel);

%% Diffusion statistics
% 12 hr traj. res.
count=0;
% f1=figure;hold on;
% f2=figure;hold on;
% cm=m_colmap('jet',360/20);

% Calculate stats as a function of time (every day for a year)
for i=1/res:1/res:size(ydisp,1)
    count=count+1;
    
    % Drift
    idx=i-((1/res)-1):i;
    ydrift(count)=mean(mean(ydisp(idx,:)/1000,'omitnan'),'omitnan');
    yhist(:,count)=histc(mean(ydisp(idx,:)/1000,'omitnan'),-300:10:400);
    
    xdrift(count)=mean(mean(xdisp(idx,:)/1000,'omitnan'),'omitnan');
    xhist(:,count)=histc(mean(xdisp(idx,:)/1000,'omitnan'),-300:10:400);
    
    % Dispersion
    qx(count)=var(mean(xdisp(idx,:)/1000,'omitnan'),'omitnan');
    qy(count)=var(mean(ydisp(idx,:)/1000,'omitnan'),'omitnan');
    
    % kurtosis
    kux(count)=kurtosis(mean(xdisp(idx,:),'omitnan')');
    kuy(count)=kurtosis(mean(ydisp(idx,:),'omitnan'));
    
%     if rem(count,20)==0
%         x=mean(xdisp(idx,:),'omitnan')';
%         x0 = x - nanmean(x);  
%         qh=histc(x0/1000,-50:5:50);
%         
%         figure(f1);
%         plot(-50:5:50,qh,'color',cm(count/20,:));
%         
%         x=mean(ydisp(idx,:),'omitnan')';
%         x0 = x - nanmean(x);  
%         qh=histc(x0/1000,-100:5:200);
%         
%         figure(f2);
%         plot(-100:5:200,qh,'color',cm(count/20,:));
%     end

    if count~=365
        mxvel(count,1)=mean(mean(xvel(idx,:)/3600,1,'omitnan').*...
            mean(xdisp(idx,:),1,'omitnan'),'omitnan');
        mxvel(count,2)=std(mean(xvel(idx,:)/3600,1,'omitnan').*...
            mean(xdisp(idx,:),1,'omitnan'),'omitnan');
        
        myvel(count,1)=mean(mean(yvel(idx,:)/3600,1,'omitnan').*...
            mean(ydisp(idx,:),1,'omitnan'),'omitnan');
        myvel(count,2)=std(mean(yvel(idx,:)/3600,1,'omitnan').*...
            mean(ydisp(idx,:),1,'omitnan'),'omitnan');
        
        mxxvel(count)=mean(mean(xvel(idx,:)/3600,1,'omitnan'),'omitnan');
        myyvel(count)=mean(mean(yvel(idx,:)/3600,1,'omitnan'),'omitnan');

    end
end

% % Diffusivity calcs- the time derivative of the second moment (dispersion)
% % change velocity units into m
% tmpx=(qx.*1000)/86400; % m2/day
% tmpy=(qy.*1000)/86400; % m/day 
% k_xx=1/2*(diff(tmpx).^2);
% k_yy=1/2*(diff(tmpy).^2);

k_xx=1/2*(diff(qx*1e6/86400)); % m^2/s
k_yy=1/2*(diff(qy*1e6/86400));

% k_xx=k_xx*1000/86400; % m^2/s
% k_yy=k_yy*1000/86400; % m^2/s

k_x=[cumtrapz(1:364,mxvel(:,1)) cumtrapz(1:364,mxvel(:,1)-mxvel(:,2)) ...
    cumtrapz(1:364,mxvel(:,1)+mxvel(:,2))];
k_y=[cumtrapz(1:364,myvel(:,1)) cumtrapz(1:364,myvel(:,1)-myvel(:,2)) ...
    cumtrapz(1:364,myvel(:,1)+myvel(:,2))];

k_xxx=var(mxxvel).*autocorr(mxvel(:,1)-mxvel(:,2));
k_yyy=var(myyvel).*autocorr(myvel(:,1)-myvel(:,2));

Tl_X=trapz(autocorr(mxvel(:,1)-mxvel(:,2)));
Tl_Y=trapz(autocorr(myvel(:,1)-myvel(:,2)));

% Tl_x=trapz(mxvel(:,1)/max(mxvel(:,1)));
% Tl_y=trapz(myvel(:,1)/max(myvel(:,1)));

% spreading stats
across_spread=mean(k_xx(1:20));
along_spread=mean(k_yy(1:100));
mx_spread=max(k_xx);
my_spread=max(k_yy);

%% DO axis stuff before plotting
fname='/ocean/rich/more/mmapbase/noaa_bc3/british_columbia_3_msl_2013.nc';
lat_lim=[48 50.3];
lon_lim=[-125.99298287 -122.22345235];
lat=ncread(fname,'lat');
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1)-.2 & lon<=lon_lim(2)+.2;
ilat=lat>=lat_lim(1)-.2 & lat<=lat_lim(2)+.2;
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);
Z(isnan(Z))=0;

PT=load('PT_AS_dist.mat');

% make axis and do rotation
lon_lim=[-127 -121];
lat_lim=[46.5 52];
m_proj('equidistant','lon',lon_lim,'lat',lat_lim);   % Projection
yl=47.85:0.01:52;
xl=ones(1,length(yl))*-122.9798;
vl=[xl;yl];
yl_cent=48.8176;xl_cent=-122.9798;
lcenter=repmat([xl_cent; yl_cent],1,length(xl));
R=[cosd(60) -sind(60);sind(60) cosd(60)];
sl=vl-lcenter;lso=R*sl;AS_axis=lso+lcenter;

% R=[-cosd(90) sind(90);-sind(90) -cosd(90)];
% sl=AS_axis-lcenter;lso=R*sl;AcS_axis=lso+lcenter;

xl=-125.4:0.01:-122.4;
yl=ones(1,length(xl))*yl_cent;
hl=[xl;yl];
lcenter=repmat([xl_cent; yl_cent],1,length(hl));
sl=hl-lcenter;lso=R*sl;AcS_axis=lso+lcenter;

AS_axis_dist=[];
AcS_axis_dist=[];
for i=1:length(AS_axis)
    AS_axis_dist(i)=gsw_distance([AS_axis(1,i) -122.9798],...
        [AS_axis(2,i) 48.8176])/1000;
    i
end

for i=1:length(AcS_axis)
    AcS_axis_dist(i)=gsw_distance([AcS_axis(1,i) -122.9798],...
        [AcS_axis(2,i) 48.8176])/1000;
    i
end

[~,idx]=min(AS_axis_dist); AS_axis_dist(1:idx)=AS_axis_dist(1:idx)*-1;
[~,idx]=min(AcS_axis_dist); AcS_axis_dist(1:idx)=AcS_axis_dist(1:idx)*-1;

%% Plotting
ylm(1,:)=[min([xdisp(:);ydisp(:)])/1000 max([xdisp(:);ydisp(:)])/1000];
ylm(2,:)=[min([qx qy]) max([qx qy])];
ylm(3,:)=[min([kux kuy]) max([kux kuy])];
ylm(4,:)=[min([k_xx';k_yy']) max([k_xx';k_yy'])];

figure('units','centimeters','outerposition',[0 0 14 15],'color','w');
ax(1)=subplot(4,3,2);
hold on
cm=viridis;cm(1,:)=[1 1 1];
colormap(ax(1),cm);
imagesc(1:length(xdrift),-300:10:400,log(xhist));
cl=caxis;
plot(1:length(xdrift),xdrift,'w','linewidth',1.75);
ylabel('Mean Drift (km)');
title('Across Strait');
axis tight; grid on; ylim(ylm(1,:));
text(0.8,0.15,'b)','units','normalized','fontweight','bold','fontsize',8,...
    'color','k');
xlim([0 200]);
xticks([0:50:200]);xticklabel_rotate([],45);


ax(2)=subplot(4,3,3);
hold on
colormap(ax(2),cm);
caxes=cl;
imagesc(1:length(ydrift),-300:10:400,log(yhist));
plot(1:length(xdrift),ydrift,'w','linewidth',1.75);
% ylabel('Mean Drift (km)');
axis tight; grid on;ylim(ylm(1,:));
title('Along Strait');
text(0.8,0.15,'c)','units','normalized','fontweight','bold','fontsize',8,...
    'color','k');
xlim([0 200]);xlim([0 200]); xticks([0:50:200]);xticklabel_rotate([],45);

c=colorbar('eastoutside');
c.Position = c.Position + [0.1 0 0.01 0];
c.Ticks=log([1 10 100 1000 5000]);
c.TickLabels=[1 10 100 1000 5000];
set(c,'tickdir','out');
title(c,'Particles','fontsize',7,'horizontalalignment','left');
% plot(1:length(xdrift),ydrift/1000,'r--');
% plot(1:length(along_drift),across_drift/1000,'b');
% plot(1:length(across_drift),along_drift/1000,'b--');

subplot(4,3,5)
plot(1:length(qx),qx,'r','linewidth',1);
ylabel('Dispersion (km^2)');
axis tight; grid on;ylim(ylm(2,:));
text(0.8,0.1,'d)','units','normalized','fontweight','bold','fontsize',8);
xlim([0 200]);xlim([0 200]); xticks([0:50:200]);xticklabel_rotate([],45);

subplot(4,3,6)
plot(1:length(qx),qy,'r','linewidth',1);
axis tight; grid on;ylim(ylm(2,:));
text(0.8,0.1,'e)','units','normalized','fontweight','bold','fontsize',8);
xlim([0 200]);xlim([0 200]); xticks([0:50:200]);xticklabel_rotate([],45);

subplot(4,3,8)
plot(1:length(kux),kux,'r','linewidth',1);
ylabel('Kurtosis');
axis tight; grid on;ylim(ylm(3,:));
text(0.8,0.1,'f)','units','normalized','fontweight','bold','fontsize',8);
set(gca,'ytick',3:6);
xlim([0 200]);xlim([0 200]); xticks([0:50:200]);xticklabel_rotate([],45);

subplot(4,3,9)
plot(1:length(kux),kuy,'r','linewidth',1);
axis tight; grid on;ylim(ylm(3,:));
text(0.8,0.1,'g)','units','normalized','fontweight','bold','fontsize',8);
set(gca,'ytick',3:6);
xlim([0 200]);xlim([0 200]); xticks([0:50:200]);xticklabel_rotate([],45);

subplot(4,3,11)
plot(1:length(k_xx),k_xx,'r','linewidth',1);
ylabel('Diffusivity (m^2 s^{-1})');
axis tight; grid on;ylim(ylm(4,:));
text(0.8,0.1,'h)','units','normalized','fontweight','bold','fontsize',8);
xlim([0 200]);xlim([0 200]); xticks([0:50:200]);xticklabel_rotate([],45);

subplot(4,3,12)
plot(1:length(k_xx),smooth(k_yy,4),'r','linewidth',1);
axis tight; grid on;ylim(ylm(4,:));
text(0.8,0.1,'i)','units','normalized','fontweight','bold','fontsize',8);
xlabel('Days since release');
% plot(1:length(k_x),k_y(:,1),'r--');
% plot(1:length(k_x),k_ac(:,1),'b');
% plot(1:length(k_x),k_al(:,1),'b--');
xlim([0 200]);xlim([0 200]); xticks([0:50:200]);xticklabel_rotate([],45);

lon_lim=[-125.5 -122.522345235];
lat_lim=[50.1 48.5];
m_proj('oblique','lon',lon_lim,'lat',lat_lim,'dir','vert','aspect',0.25)
ax2=axes('Position',[0.1 0.11 .2 0.8]);
% caxis([0 400]);
cm=flipud(m_colmap('blues',50)); cm=cm(end-40:end,:);
[CS,h]=m_contourf(lon(ilon),lat(ilat),Z'*-1,0:25:400,'edgecolor','none');
colormap(ax2,cm);
hold on
m_gshhs_f('patch',rgb('light grey'),'edgecolor','none');
m_northarrow(-125.32,49.75,0.15,'type',2);
m_plot(PT.mn_traj(1:24*100,1),PT.mn_traj(1:24*100,2),'color',[1 1 1],...
    'linewidth',1.75);
m_grid('linestyle','none','tickdir','out','xaxisloc','bottom','yaxisloc','left','fontsize',8,...
    'xtick',-126:-122,'ytick',48:51);
% Plot representative stations

m_line(AcS_axis(1,200:end-25),AcS_axis(2,200:end-25),'color',rgb_x('light red'),'linewidth',2);
m_annotation('arrow',[AcS_axis(1,end-25),AcS_axis(1,end-24)],...
    [AcS_axis(2,end-25),AcS_axis(2,end-24)],'color',rgb_x('light red'),'linewidth',1);

m_line(AS_axis(1,50:end-50),AS_axis(2,50:end-50),'color',rgb_x('light red'),'linewidth',2);
m_annotation('arrow',[AS_axis(1,end-50),AS_axis(1,end-49)],...
    [AS_axis(2,end-50),AS_axis(2,end-49)],'color',rgb_x('light red'),'linewidth',1);

m_scatter(-122.9798,48.8176,75,'w','filled');
m_scatter(-122.9798,48.8176,50,'x','k','linewidth',1);
text(0.75,0.05,'a)','units','normalized','fontweight','bold','fontsize',8);
c2=colorbar('southoutside');
c2.Position = c2.Position + [-0.025 -0.125 0.05 0];
c2.Ticks=0:200:400;c2.TickLabels={'0 m';'200 m';'400 m'};
plot(station_bound,'FaceColor','red','FaceAlpha',0.1,'edgecolor','none')
set(findall(gcf,'-property','FontSize'),'FontSize',8);

count=-1;
for i=-75:25:250
    [~,idx]=min(abs(AS_axis_dist-i));
    m_scatter(AS_axis(1,idx),AS_axis(2,idx),20,'filled',...
        'markerfacecolor','w','markeredgecolor','r');
    if rem(count,2)==0 && AS_axis(1,idx)>lon_lim(1) && ...
            AS_axis(1,idx)<lon_lim(2) && ...
            AS_axis(2,idx)>lat_lim(2) && ...
            AS_axis(2,idx)<lat_lim(1) && ...
            i~=0
        
        m_text(AS_axis(1,idx)+0.01,AS_axis(2,idx)+0.01,...
            sprintf('%3.0f km',i),'FontSize',6,'color','r');
    end
    count=count+1;
end

for i=-75:25:max(AcS_axis_dist)
    [~,idx]=min(abs(AcS_axis_dist-i));
    m_scatter(AcS_axis(1,idx),AcS_axis(2,idx),20,'filled',...
        'markerfacecolor','w','markeredgecolor','r');
    if i==-25 || i==25
        m_text(AcS_axis(1,idx)+0.05,AcS_axis(2,idx)-0.05,...
            sprintf('%3.0f km',i),'FontSize',6,'color','r');
        %         m_annotation('textbox',[AcS_axis(1,idx)+0.05 AcS_axis(2,idx)-0.05 0.1 0.1],...
        %             'String',sprintf('%3.0f km',i),'FontSize',6,'color','k','FitBoxToText','on');
    end
end

set(findall(gcf,'-property','FontSize'),'FontSize',8)

% subplot(2,3,3)
% plot(1:length(kux),kux,'k');
% hold on
% plot(1:length(kux),kuy,'k--');
% plot(1:length(kux),kacross,'r');
% plot(1:length(kux),kalong,'r--');
% ylabel('Kurtosis');
% axis tight; grid on;
% % Initially, the kurtosis values are high as particles are drawn out in
% % filaments. These values decrease to a minima at time tm after, indicating
% % the initial dominance of diffusion on short time scales. The timing of
% % this minima is not uniform for all directions.  Then, the kurtosis
% % increases to a maxima as coherent advective structures distribute
% % particles. If the motion was mainly diffusion, the kurtosis value would stay low.
%
% % subplot(2,3,2)
% % plot(res:res:365-res,absxdiff,'k');
% % hold on
% % plot(res:res:365-res,absydiff,'r');
% % legend('Across-strait','Along-strait','location','best');
% % % ylabel(a
% % axis tight; grid on;


%%
export_fig /ocean/sstevens/IW_project/figures/paper/vec/PT_stats_V3.png -m3
export_fig /ocean/sstevens/IW_project/figures/paper/vec/PT_stats_V3.pdf -dpdf

% %%
% all_pass.lon=[];
% all_pass.lat=[];
% all_pass.age=[];
% all_pass.z=[];
% 
% % compile pass info
% for i=1:length(pass_info)
%     
%     nanmsk=isnan(pass_info(i).lon);
%     pass_info(i).lon(nanmsk)=[];
%     pass_info(i).lat(nanmsk)=[];
%     pass_info(i).age(nanmsk)=[];
%     pass_info(i).z(nanmsk)=[];
%     
%     all_pass.lon=[all_pass.lon;pass_info(i).lon];
%     all_pass.lat=[all_pass.lat;pass_info(i).lat];
%     all_pass.age=[all_pass.age;pass_info(i).age];
%     all_pass.z=[all_pass.z;pass_info(i).z];
% end
% 
% % Sort by age
% all_pass.age(all_pass.age==-1)=NaN;
% [all_pass.age,idx]=sort(all_pass.age);
% all_pass.lon=all_pass.lon(idx);
% all_pass.lat=all_pass.lat(idx);
% all_pass.z=all_pass.z(idx);
% 
% % find particles that pass through EBC
% all_pass.lon(all_pass.lon==-1)=NaN;
% 
% % EBC core
% EBC_pass=all_pass.lon>-123.3088;
% 
% % % Whole EBC
% % EBC_pass=all_pass.lon>-123.3233;
% 
% gate_pass_CP=[];
% for i=1:365
%     gate_pass_CP=[gate_pass_CP;sum(...
%         all_pass.age>i-1 & all_pass.age<i)];
% end
% gate_pass_CP(:,2)=(cumsum(gate_pass_CP(:,1))/length(all_pass.age))*100;
% 
% EBC_CP=[];
% for i=1:365
%     EBC_CP=[EBC_CP;sum(...
%         all_pass.age(EBC_pass)>i-1 & all_pass.age(EBC_pass)<i)];
% end
% EBC_CP(:,2)=(cumsum(EBC_CP(:,1))/length(all_pass.age))*100;
% 
% figure('units','centimeters','outerposition',[0 0 11.5 11.5],'color','w');
% histogram(all_pass.age,0:2.5:150,'facecolor',[1 1 1]*0.8,'edgecolor','k');
% hold on
% ax=histogram(all_pass.age(EBC_pass),0:2.5:150,'facecolor',rgb_x('light red')...
%     ,'edgecolor','k');
% xlim([0 75]);
% % ytick(0:1000:5000);
% ylabel('Number of particles N','fontweight','bold');
% xlabel('Transit time (days)','fontweight','bold');
% yyaxis right
% plot(1:365,gate_pass_CP(:,2),'k');
% hold on
% plot(1:365,EBC_CP(:,2),'color',rgb_x('light red'));
% ylabel('Cumulative % N','fontweight','bold');
% set(gca,'ycolor','k');
% set(findall(gcf,'-property','FontSize'),'FontSize',7)
% 
% % % set(gca,'XScale','log')
% %
% % figure('units','centimeters','outerposition',[0 0 11.5 11.5],'color','w');
% % histogram(all_pass.age,1:0.25:10);
% % hold on
% % histogram(all_pass.age(EBC_pass),1:0.25:10);
% % histogram(all_pass.age(~EBC_pass),1:0.25:10);
% 
% export_fig '/ocean/sstevens/IW_project/figures/paper/PT_pass.eps' -depsc
% export_fig '/ocean/sstevens/IW_project/figures/poster/PT_passPOSTER.png' -png -m3 -transparent


