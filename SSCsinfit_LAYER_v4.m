%% Do fitting for all depth levels in SSC 
% surface
clc
clear
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
load /ocean/sstevens/IW_project/data/model_IW_outputT_surf.mat
close
load /ocean/sstevens/IW_project/data/no_inlet_msk.mat
load /ocean/sstevens/IW_project/data/thalweg.mat
lnes=lines;
grey=rgb('light grey');
model_temp.lon=double(model_temp.lon);model_temp.lat=double(model_temp.lat);
obs=load('sinfit_temp_data.mat','station_lon','station_lat');

load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;

inside_idx=inShape(station_bound,model_temp.lon,model_temp.lat);
grid_lon=model_temp.lon; grid_lon(~inside_idx)=NaN;
grid_lat=model_temp.lat; grid_lat(~inside_idx)=NaN;

%% Create sinusoid fits to find coldest yearday
cold_day=NaN(size(grid_lon));
% fitted_harms=NaN([size(grid_lon) length(0.01:0.01:365)]);
amp=NaN(size(grid_lon));

count=0;
tmp_day=day(datetime(datestr(model_temp.mtime)),'dayofyear');
hh = waitbar(0,'Fitting sinusoids...');
harm_vec=0.01:0.01:365;
for i=1:size(grid_lon,1)
    for ii=1:size(grid_lon,2)
        count=count+1;
        if ~isnan(grid_lon(i,ii)) && sum(model_temp.t(i,ii,:))>0
            
            % Pinpoint closest model_temp points to stations
            tmp_temp=squeeze(model_temp.t(i,ii,:));
            
            % fit sinusoid
            [amp(i,ii),phase,frac,offset,yy]=fit_harmonics(tmp_temp,tmp_day',...
                1,365,0.01);
            
            count2=0;
            for tt=0.01:0.01:365
                count2=count2+1;
                fitted_harm(count2)=offset+amp(i,ii)*cos(2*pi*tt/365 + phase);
            end
            
            [~,idx]=min(fitted_harm);
            cold_day(i,ii)=harm_vec(idx);
%             fitted_harms(i,ii,:)=fitted_harm;
        end
            waitbar(count/numel(grid_lon),hh);
    end
end

% Save data for plotting elsewhere
clearvars -except cold_day fitted_harms amp grid_lon grid_lat
save('/ocean/sstevens/IW_project/data/sinfit_SSC_modelres_surf.mat');

%% core
clear
load /ocean/sstevens/IW_project/data/model_IW_outputT_core.mat
close
load /ocean/sstevens/IW_project/data/no_inlet_msk.mat
load /ocean/sstevens/IW_project/data/thalweg.mat
lnes=lines;
grey=rgb('light grey');
model_temp.lon=double(model_temp.lon);model_temp.lat=double(model_temp.lat);
obs=load('sinfit_temp_data.mat','station_lon','station_lat');

load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;

inside_idx=inShape(station_bound,model_temp.lon,model_temp.lat);
grid_lon=model_temp.lon; grid_lon(~inside_idx)=NaN;
grid_lat=model_temp.lat; grid_lat(~inside_idx)=NaN;

%% Create sinusoid fits to find coldest yearday
cold_day=NaN(size(grid_lon));
fitted_harms=NaN([size(grid_lon) 365]);
amp=NaN(size(grid_lon));

count=0;
tmp_day=day(datetime(datestr(model_temp.mtime)),'dayofyear');
hh = waitbar(0,'Fitting sinusoids...');
for i=1:size(grid_lon,1)
    for ii=1:size(grid_lon,2)
        count=count+1;
        if ~isnan(grid_lon(i,ii)) && sum(model_temp.t(i,ii,:))>0
            
            % Pinpoint closest model_temp points to stations
            tmp_temp=squeeze(model_temp.t(i,ii,:));
            
            % fit sinusoid
            [amp(i,ii),phase,frac,offset,yy]=fit_harmonics(tmp_temp,tmp_day',...
                1,365,0.01);
            
            for tt=1:365
                fitted_harm(tt)=offset+amp(i,ii)*cos(2*pi*tt/365 + phase);
            end
            
            [~,cold_day(i,ii)]=min(fitted_harm);
            fitted_harms(i,ii,:)=fitted_harm;
        end
            waitbar(count/numel(grid_lon),hh);
    end
end

% Save data for plotting elsewhere
clearvars -except cold_day fitted_harms amp grid_lon grid_lat
save('/ocean/sstevens/IW_project/data/sinfit_SSC_modelres_core.mat');

%% deep
clear
load /ocean/sstevens/IW_project/data/model_IW_outputT_deep.mat
close
load /ocean/sstevens/IW_project/data/no_inlet_msk.mat
load /ocean/sstevens/IW_project/data/thalweg.mat
lnes=lines;
grey=rgb('light grey');
model_temp.lon=double(model_temp.lon);model_temp.lat=double(model_temp.lat);
obs=load('sinfit_temp_data.mat','station_lon','station_lat');

load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;

inside_idx=inShape(station_bound,model_temp.lon,model_temp.lat);
grid_lon=model_temp.lon; grid_lon(~inside_idx)=NaN;
grid_lat=model_temp.lat; grid_lat(~inside_idx)=NaN;

%% Create sinusoid fits to find coldest yearday
cold_day=NaN(size(grid_lon));
fitted_harms=NaN([size(grid_lon) 365]);
amp=NaN(size(grid_lon));

count=0;
tmp_day=day(datetime(datestr(model_temp.mtime)),'dayofyear');
hh = waitbar(0,'Fitting sinusoids...');
for i=1:size(grid_lon,1)
    for ii=1:size(grid_lon,2)
        count=count+1;
        if ~isnan(grid_lon(i,ii)) && sum(model_temp.t(i,ii,:))>0
            
            % Pinpoint closest model_temp points to stations
            tmp_temp=squeeze(model_temp.t(i,ii,:));
            
            % fit sinusoid
            [amp(i,ii),phase,frac,offset,yy]=fit_harmonics(tmp_temp,tmp_day',...
                1,365,0.01);
            
            for tt=1:365
                fitted_harm(tt)=offset+amp(i,ii)*cos(2*pi*tt/365 + phase);
            end
            
            [~,cold_day(i,ii)]=min(fitted_harm);
            fitted_harms(i,ii,:)=fitted_harm;
        end
            waitbar(count/numel(grid_lon),hh);
    end
end

% Save data for plotting elsewhere
clearvars -except cold_day fitted_harms amp grid_lon grid_lat
save('/ocean/sstevens/IW_project/data/sinfit_SSC_modelres_deep.mat');
