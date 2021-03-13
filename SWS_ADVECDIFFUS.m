% Advection/Diffusion equation solutions
addpath(genpath('/ocean/sstevens/'));
addpath(genpath('/ocean/rich/home/matlab/m_map/'));
%%
clear

x=[0:.01:100];

K=2*(sqrt(1+i*(2*x).^2)-1)./(2*x).^2;

h=semilogx(x,real(K),x,imag(K));
set(h(1),'linewi',2);
legend('k_R','k_I');
xlabel('X_D/X_A');
set(gca,'tickdir','out','fontsize',16,'box','off');

% %print -depsc Advecdiffus
% print -dpng Advecdiffus

%%
t=linspace(0,2*pi,80);
x=0:.01:1;
[xx,tt]=meshgrid(x,t);

XA=1/(2*pi)*2;  % ADVECTIVE LENGTH SCALE
                % "2" means twice length in a year = "Strait of Georgia"

Gam=25;
Gam=.25;
Gam=2.5;   % 2*XD/XA
           % = 2.5 for "Strait of Georgia"

% Nondimensional soln.

K=2/XA*(sqrt(1+i*Gam.^2)-1)./Gam.^2;

 
S=-real(exp(i*tt-K*xx));

clf;set(gcf,'defaultaxestickdir','out','defaultaxestickdirmode','manual');
subplot(2,1,1);
hh=plot(x,S(1,:),x,0*x,'k',x,exp(-real(K)*x),'--',x,-exp(-real(K)*x),'--');
set(hh(1),'linewi',2);
title('Envelope and solution at t=0');
xlabel('X');

subplot(2,1,2);
%imagesc(x,t/(2*pi)*365,S);set(gca,'ydir','normal');
contourf(x,t/(2*pi)*365,S,[-1:.1:1]);
colormap(jet);
xlabel('X (normalized)');ylabel('T (days)');
line([0 1],[0 1/XA]/(2*pi)*365,'color','w','linewi',4);  % ADVECTION SPEED LINE
line([0 1],[0 imag(K)]/(2*pi)*365,'color','w','linewi',2,'linest','--');  % Effective line
title(['Time/Range plot - X_D/X_A=' num2str(Gam/2)]);

% print('-dpng',['Advecdiffus2_' num2str(Gam/2) '.png']);


%% Solutions for delta function

%% Do this again...For real!!!!

w=2*pi/(365); %rad/day
X=200;  % km
% U=X/(.5*365);  % km/day
U=2e-5*(60*60*24); % 1.5 cm/s (in km/day) as per measurements
% % U=1e-5*(60*60*24); % 1 cm/s
% U=0.5e-5*(60*60*24); % 0.5 cm/s

% gam = 2A/U
Gam=[20 1e2 1e4];

figure('units','centimeters','outerposition',[0 0 14 15],'color','w');
% orient portrait;wysiwyg;
set(gcf,'defaultaxestickdir','out','defaultaxestickdirmode','manual');

count=0;
lab={'d)' 'a)' 'e)' 'b)' 'f)' 'c)'};
for kk=1:3
count=count+1;
k=(sqrt(1+2*i*w/U*Gam(kk))-1)/Gam(kk);

t=linspace(0,365,80);
x=linspace(0,X,200);
[xx,tt]=meshgrid(x,t);

S=-real(exp(-k*xx+i*w*tt));

ax(count)=subplot(3,3,3+kk);
xlim([-1.1 1.1]);
hold on
hh=plot(S(1,:),x,0*x,x,'k:',exp(-real(k)*x),x,'--',-exp(-real(k)*x),x,'--');
set(hh(1),'linewi',1.5,'color','k');
set(hh(3),'color',rgb_x('light red'));
set(hh(4),'color',rgb_x('light red'));
patch([exp(-real(k)*x) fliplr(-exp(-real(k)*x))],[x fliplr(x)],rgb_x('light red'),...
    'facealpha',0.2,'linestyle','none');
text(0.05,1.1,lab{count},'units','normalized');

count=count+1;
ax(count)=subplot(3,3,kk);
%imagesc(x,t/(2*pi)*365,S);set(gca,'ydir','normal');
[CS,CH]=contourf(t,x,fliplr(rot90(S,-1)),-1:0.2:1,'linestyle','none');%,[-1:.1:1]);
caxis([-1 1]);
colormap(cmocean('balance','pivot',0));
% xlabel('X  (km)');ylabel('Time (days)');
line([0 X/U],[0 X],'color','k','linewi',3);  % ADVECTION SPEED LINE
line([0 X/U],[0 X],'color','w','linewi',1.5);  % ADVECTION SPEED LINE
line([0 imag(k)*X/w],[0 X],'color','w','linewi',2,'linest','--');  % Effective line
% title(['Time/Range plot - 2A/U=' num2str(Gam(kk))]);
title(['\gamma=' num2str(Gam(kk)) ' km']);
set(gca,'xtick',0:100:300);
text(0.05,1.1,lab{count},'units','normalized');

end

axes(ax(2))
ylabel('Distance (km)','fontsize',8);
axes(ax(1))
ylabel('Distance (km)','fontsize',8);
axes(ax(3))
xlabel('Temperature (^\circC)','fontsize',8);
axes(ax(4))
xlabel('Time (days)','fontsize',8);

axes(ax(6))
[ax,h]=m_contfbar([0.25 0.75],-.2,CS,CH,'endpiece','no','axfrac',.075,...
    'fontsize',6,'linest','none');

Gam=logspace(-3,6,100);
k=(sqrt(1+2*i*w/U.*Gam)-1)./Gam;
lnes=lines;

subplot(3,3,7:9);
h=semilogx(Gam,real(k)*X,Gam,imag(k)*X);
xlim(Gam([1 end]));
set(h(1),'linewi',1,'color',lnes(1,:));
set(h(2),'color',lnes(2,:));
grid on
hold on
patch([56 56 75 75],[0 max(imag(k)*X) max(imag(k)*X) 0],...
    rgb_x('light blue'),...
    'facealpha',0.4,'linestyle','none');
line([20 20],[0 max(imag(k)*X)],'color','k','linewi',0.5); 
line([100 100],[0 max(imag(k)*X)],'color','k','linewi',0.5);
line([1e4 1e4],[0 max(imag(k)*X)],'color','k','linewi',0.5);
xlabel('\gamma (km)');
ylabel('\it{k} \rm(km ^{-1})');
set(gca,'tickdir','out','fontsize',8,'box','off');
legend(h,'\it{k_r}','\it{k_i}');
text(0.02,1.1,'g)','units','normalized');

set(findall(gcf,'-property','FontSize'),'FontSize',8);

%% 
export_fig /ocean/sstevens/IW_project/figures/paper/vec/advecdiffus_V3.pdf -dpdf -nocrop

%%


figure(2);
clf;set(gcf,'defaultaxestickdir','out','defaultaxestickdirmode','manual');

Gam=logspace(-3,6,100);
k=(sqrt(1+2*i*w/U.*Gam)-1)./Gam;

h=semilogx(Gam,real(k)*X,Gam,imag(k)*X);
xlim(Gam([1 end]));
set(h(1),'linewi',1);
legend('k_RX','k_IX');
xlabel('2A/U (km)');
set(gca,'tickdir','out','fontsize',16,'box','off');

%%
figure(1);print -dpng fosc
figure(2);print -dpng fkvals

%%

Avals=[10 100 1000 10000]; % m^2/s

t=linspace(0,3,180); % years
x=linspace(0,200,200); % km
U=200/.5; % km/year
[xx,tt]=meshgrid(x,t);

figure(1);clf;orient landscape;wysiwyg;
for kk=1:length(Avals)

    subplot(2,4,kk);  
    A=Avals(kk)*(1e-6*(365*86400)); % Convert from m^2/s to km^2/year

    c=1./sqrt(4*pi*A*tt).*exp( -(xx-U*tt).^2./(4*A*tt) );
    c(1,:)=0;
    
    contourf(x,t,(c),20);
    line(x,x./U,'color','w','linewi',2);
    line(x,x./U,'color','r','linest','--');
    caxis([0 max(c(:,1))]);
    ylim([0 1]);
    xlabel('x/km');
    if kk==1, ylabel('time/years'); end
    title(sprintf('A=%g m^2/s',Avals(kk)));
    
    subplot(2,4,4+kk);
    plot(t,c(:,end));
    if kk==1, yy=ylim; end;
    line(x(end)./[U U],yy,'color','r','linest','--');
    mn=sum(t*c(:,end))/sum(c(:,end));
    line([mn mn],yy,'color',[0 .5 0],'linewi',3);
    text(mn,yy(2)/2,sprintf('Mean %.1f',mn));
    xlabel('time/years');
    ylabel('c(200 km,t)');
    
end
colormap(m_colmap('jet'));

% print -dpng Hdiffus
    


