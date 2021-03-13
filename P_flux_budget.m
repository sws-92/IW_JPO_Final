addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));

%% Load bathy data
clear
dir_list=dir('/ocean/sstevens/ariane_old/results_BP/2016*365d');

tic

res=1; % Daily res
dep=single(NaN(length(res:res:365),20352));

load PT_AS_inlets.mat
station_bound=alphaShape(AS.lon,AS.lat);
station_bound.Alpha=station_bound.Alpha*2;

load HS_AS.mat
station_bound2=alphaShape(l,ll);

count2=0;
for i=1:length(dir_list)
    
    % Load data
    path = [dir_list(i).folder '/'  dir_list(i).name '/'];
    filestr=[path 'TrajTBL'];
    load(filestr);
    disp(i)
    
    for j=1:length(TrajTBL)
        mn_traj=NaN(length(res:res:365),3);
        count=0;
        for ii=res:res:365
            count=count+1;
            mask=(TrajTBL(j).Age_days>ii-res & TrajTBL(j).Age_days<ii);
            s(count)=sum(mask);
            mn_traj(count,1)=mean(TrajTBL(j).lon(mask));
            mn_traj(count,2)=mean(TrajTBL(j).lat(mask));
            mn_traj(count,3)=mean(TrajTBL(j).z_m(mask));
        end
        
        msk=isnan(mn_traj(:,1));
        if sum(msk)
            mn_traj(msk,:)=0;
        end
        
%         msk=~inShape(station_bound,mn_traj(:,1),mn_traj(:,2));
%         mn_traj(msk,:)=NaN;
        
        msk=~inShape(station_bound,mn_traj(:,1),mn_traj(:,2)) |...
            inShape(station_bound2,mn_traj(:,1),mn_traj(:,2));
        mn_traj(msk,:)=NaN;

        count2=count2+1;
        dep(:,count2)=mn_traj(:,3);
    end
end

toc
clearvars TrajTBL 

save('/ocean/sstevens/ariane_old/results_BP/mn_dep_excl_HS.mat');

%% Calculate flux through boundaries
clear
load /ocean/sstevens/ariane_old/results_BP/mn_dep_excl_HS.mat

dep=dep*-1;

surface=dep-50;
deep=dep-150;
rMask=dep<50 & dep>150;rIW=dep;rIW(rMask)=NaN;
pRemain=sum(~isnan(rIW),2);

sFlip=-1*diff(sign(surface))/2;
sFlux=nansum(sFlip,2);
sCum=cumtrapz(sFlux);

tmp=sFlip;msk=tmp==-1;tmp(msk)=NaN;
sIn=nansum(tmp,2);
tmp=sFlip;tmp(~msk)=NaN;
sOut=abs(nansum(tmp,2));

dFlip=diff(sign(deep))/2;
dFlux=nansum(dFlip,2);
dCum=cumtrapz(dFlux);

tmp=dFlip;msk=tmp==-1;tmp(msk)=NaN;
dIn=nansum(tmp,2);
tmp=dFlip;tmp(~msk)=NaN;
dOut=abs(nansum(tmp,2));

inFlux=dIn+sIn;
inRat=inFlux./pRemain(1:end-1);

all_flux=(sFlux+dFlux);
rat=all_flux./pRemain(2:end);
weekFlux=trapz(rat(1:7));
rCum=trapz(rat(1:150));
% ef=interp1(cumtrapz(rat),1:364,0.63);
S2D=trapz(sFlux)/trapz(dFlux);
S2Drun=nanmean(sFlux./dFlux);

Ft=interp1(cumtrapz(rat),1:364,0.63);