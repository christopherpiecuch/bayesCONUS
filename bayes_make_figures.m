%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function bayes_make_figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original code written by CGP 09 August 2025
% Updated code written by CGP 19 October 2025
% This code makes the figures in Piecuch (2025)
% NB factors of 1e3 are to convert from m to mm and
% factors of 1e2 are to convert from m to cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayes_make_figures(experimentName,runnum)

% load all gridded solution files
a=[]; b=[]; bg=[]; l=[]; y=[]; yg=[]; bg=[]; gg=[];
for kk=1:numel(runnum), disp(num2str(kk))
    load(['bayes_model_solutions/experiment_',experimentName,'_runNum_',num2str(runnum(kk)),'_gridded.mat']);
    a=[a; A];
    l=[l; L];
    y=[y; Y];
    yg=[yg; YG];
    bg=[bg; BG];
    gg=[gg; GG];
    clear A L Y YG B BG GG G
end
A=a; L=l; Y=y; YG=yg; clear a l y yg
BG=bg; GG=gg; clear bg gg
NumIter=numel(A)/numel(LON);

% load color map and define example locations
load clr0.mat
T=1900:2024; T=T-mean(T);
exlocs=[46 50 14 15 31 27]; % example locations

% predict data from model solution (Equation 8 in Piecuch et al. 2018)
LD=permute(reshape(reshape(L,NumIter*numel(NAME),1)*ones(1,numel(T)),NumIter,numel(LON),numel(T)),[1 3 2]);
AD=permute(reshape(reshape(A,NumIter*numel(NAME),1)*T,NumIter,numel(LON),numel(T)),[1 3 2]);
YD=Y+LD+AD;

% define instantaneous quadratic and linear rates
RATEGQ=permute(reshape(2*reshape(GG,numel(GG),1)*T+...
        reshape(BG,numel(BG),1)*ones(size(T)),NumIter,numel(GLON),numel(T)),[1 3 2]);
RATEGL=permute(reshape(reshape(BG,numel(BG),1)*ones(size(T)),NumIter,numel(GLON),numel(T)),[1 3 2]);

% define quadratic and linear fits
QuadPart=permute(reshape(reshape(GG,numel(GG),1)*T.^2,NumIter,numel(GLON),numel(T)),[1 3 2]); % quadratic part of solution
LinPart=permute(reshape(reshape(BG,numel(BG),1)*T,NumIter,numel(GLON),numel(T)),[1 3 2]); % linear part of solution
QuadPartConst=squeeze(mean(QuadPart,2)); % constant offset in quadratic part (positive definite) to add as intercept to linear part
QuadPartConst=reshape(reshape(QuadPartConst,numel(QuadPartConst),1)*ones(size(T)),NumIter,numel(GLON),numel(T));
QuadPartConst=permute(QuadPartConst,[1 3 2]);
SLevQ=QuadPart+LinPart;
SLevL=QuadPartConst+LinPart;

% make prelim figure to get lon/lat bounds
figure
plot(GLON,GLAT,'k.')
plot(LON,LAT,'k.')
ax=axis;
close all

% make figure 1
fig=figure('color','white');
fig.Position(4) = fig.Position(4)*1.25;
fig.Position(3) = fig.Position(3)*1.5;

 sp1=subplot(2,7,1:4);
 worldmap([20 50],[-130 -60])
 geoshow('landareas.shp','facecolor',[.9 .9 .9],'edgecolor',[.5 .5 .5])
 hold on, box on, grid on
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow(states,'facecolor',[.9 .9 .9],'edgecolor',[.5 .5 .5])
 scatterm(LAT,LON,30,sum(~isnan(DATA')),'filled','markeredgecolor',[.5 .5 .5]), caxis([30 100]),
 colormap(sp1,parula(14))
 c=colorbar('location','eastoutside');
 cpos=get(c,'Position');
 cpos(2)=cpos(2)+0.075;
 cpos(1)=cpos(1)-0.05;
 cpos(3)=cpos(3)*0.5;
 cpos(4)=cpos(4)-0.175;
 set(c,'Position',cpos)
 title('a. Years of Data at Tide Gauges','fontweight','normal')
 set(gca,'fontsize',12)

 sp2=subplot(2,7,8:11);
 worldmap([20 50],[-130 -60])
 geoshow('landareas.shp','facecolor',[.9 .9 .9],'edgecolor',[.5 .5 .5])
 hold on, box on, grid on
 states = shaperead('usastatelo', 'UseGeoCoords', true);
 geoshow(states,'facecolor',[.9 .9 .9],'edgecolor',[.5 .5 .5])
 scatterm(GLAT,GLON,30,median(squeeze(YG(:,end,:)-YG(:,1,:))),'filled','s')
 scatterm(LAT,LON,30,median(squeeze(Y(:,end,:)-Y(:,1,:))),'filled','markeredgecolor',[.5 .5 .5]), 
 colormap(sp2,turbo(20))
 c=colorbar('location','eastoutside');
 cpos=get(c,'Position');
 cpos(2)=cpos(2)+0.075;
 cpos(1)=cpos(1)-0.05;
 cpos(3)=cpos(3)*0.5;
 cpos(4)=cpos(4)-0.175;
 set(c,'Position',cpos)
 title('b. Sea-Level Rise 1900-2024 (m)','fontweight','normal')
 set(gca,'fontsize',12)
 caxis([0 1])

 % the "ABC+kk*3/10" factor is for vertical offsetting of time series
 subplot(1,7,6:7)
 hold on, box on, grid off
 for kk=1:numel(exlocs)
  ABC=-nanmean(prctile(squeeze(YD(:,:,exlocs(kk))),50));
  c=fill([1900:2024 2024:-1:1900],ABC+kk*3/10+[prctile(squeeze(YD(:,:,exlocs(kk))),2.5) fliplr(prctile(squeeze(YD(:,:,exlocs(kk))),97.5))],[.5 .5 .5]);
  set(c,'edgecolor','none','FaceAlpha',0.5)
  plot(1900:2024,ABC+kk*3/10+prctile(squeeze(YD(:,:,exlocs(kk))),50),'k','linewidth',2)
 plot(1900:2024,ABC+kk*3/10+DATA(exlocs(kk),:),'.-','color',clr0(7,:),'markersize',5,'linewidth',2)
end
axis([1900 2025 0 2.2])
set(gca,'ytick',(1:6)*3/10-0.1,'yticklabel',[{'New York, NY'};{'Nantucket, MA'};{'Fernandina Beach, FL'};{'Fort Pulaski, GA'};{'San Francisco, CA'};{'Monterey, CA'}])
ytickangle(25)
set(gca,'fontsize',12,'fontweight','normal','xtick',1900:25:2025)
plot([2000 2000],[0.05 0.35],'k','linewidth',2)
plot([1999 2001],[0.05 0.05],'k','linewidth',2)
plot([1999 2001],[0.35 0.35],'k','linewidth',2)
text(1998,0.2,'30 cm','fontsize',12,'fontweight','normal','HorizontalAlignment','right')
legend([{'This Study (95% CI)'};{'This Study (Median)'};{'Tide Gauge'}],'location','northwest','orientation','vertical')
legend boxoff
title([{'c. RSL Changes at'};{'Example Locations'}],'fontweight','normal')
xlabel('Year','fontsize',12,'fontweight','normal')

% define cos(latitude) weighting mask
ijk=find(GLON>=-100&GLON<=-90); % define the western gulf
weight=reshape(ones(NumIter,1)*cos(pi/180*GLAT'),NumIter,numel(GLON));
weight2=weight; % this is the mask for the no-western-gulf calculation
weight2(:,ijk)=0; % mask for omitting western gulf
WEIGHT=reshape(ones(NumIter*125,1)*cos(pi/180*GLAT'),NumIter,numel(T),numel(GLON));
WEIGHT2=WEIGHT; % this is the mask for the no-western-gulf calculation
WEIGHT2(:,:,ijk)=0; % mask for omitting western gulf

% compute conus averages
YGM=sum(YG.*WEIGHT,3)./sum(WEIGHT,3);
RATEQ=sum(RATEGQ.*WEIGHT,3)./sum(WEIGHT,3);
RATEL=sum(RATEGL.*WEIGHT,3)./sum(WEIGHT,3);
ACCQ=sum(2*GG.*weight,2)./sum(weight,2);
TRENDL=sum(BG.*weight,2)./sum(weight,2);
SLQ=sum(SLevQ.*WEIGHT,3)./sum(WEIGHT,3);
SLL=sum(SLevL.*WEIGHT,3)./sum(WEIGHT,3);

% compute conus averages less the western gulf
YGM2=sum(YG.*WEIGHT2,3)./sum(WEIGHT2,3);
RATEQ2=sum(RATEGQ.*WEIGHT2,3)./sum(WEIGHT2,3);
RATEL2=sum(RATEGL.*WEIGHT2,3)./sum(WEIGHT2,3);
ACCQ2=sum(2*GG.*weight2,2)./sum(weight2,2);
TRENDL2=sum(BG.*weight2,2)./sum(weight2,2);
SLQ2=sum(SLevQ.*WEIGHT2,3)./sum(WEIGHT2,3);
SLL2=sum(SLevL.*WEIGHT2,3)./sum(WEIGHT2,3);

% compute 25-year residuals
reswin=25;
RESL=squeeze(mean(reshape(YGM-SLL,NumIter,reswin,numel(YEAR)/reswin),2));

% figure 2
fig=figure('color','white');
fig.Position(3) = fig.Position(3)*1.5;
fig.Position(4) = fig.Position(4)*1.5;
 subplot(3,3,[1 2 4 5])
 hold on, box on, grid on
 plot(1900:2024,1e2*prctile(YGM,50),'k','linewidth',2)
 c=fill([1900:2024 2024:-1:1900],1e2*[prctile(YGM,2.5) fliplr(prctile(YGM,97.5))],[.5 .5 .5]);
 set(c,'edgecolor','none','FaceAlpha',0.5)

 plot(1900:2024,1e2*prctile(SLQ,50),'color',clr0(2,:),'linewidth',2)
 plot(1900:2024,1e2*prctile(SLQ,97.5),'--','color',clr0(2,:),'linewidth',1)
 plot(1900:2024,1e2*prctile(SLQ,2.5),'--','color',clr0(2,:),'linewidth',1)
 axis([1900 2024 -30 30])
 legend([{'CONUS RSL (median)'};{'CONUS RSL (95% CI)'};{'Quadratic (median)'};{'Quadratic (95% CI)'}],'location','northwest'), legend boxoff
 axis([1900 2025 -30 30])
 ylabel('RSL (cm)','fontsize',12,'fontweight','normal')
 xlabel('Year','fontsize',12,'fontweight','normal')
 set(gca,'fontsize',12)
 title('a. CONUS RSL (cm)','fontweight','normal')
 set(gca,'fontsize',12,'fontweight','normal','xtick',1900:25:2025,'ytick',-30:10:30)
 set(gca,'fontsize',12)

 subplot(3,3,[7 8])
 b=bar(1:5,1e2*median(RESL),'facecolor',clr0(1,:)*0.5+[.5 .5 .5],'edgecolor',clr0(1,:));               
 hold on
 er = errorbar(1:5,1e2*median(RESL),1e2*prctile(RESL,2.5),1e2*prctile(RESL,97.5));    
 er.Color = clr0(1,:);                            
 er.LineStyle = 'none';  
 er.LineWidth = 1;  
 set(gca,'xtick',1:5,'xticklabel',[{'1900-1924'};{'1925-1949'};{'1950-1974'};{'1975-1999'};{'2000-2024'}])
 xtickangle(45)
 axis([0 6 -2 2])
 ylabel('25-Year Residual (cm)','fontweight','normal')
 title([{'b. CONUS RSL 25-Year-Average Residuals From Historical Trend (cm)'}],'fontweight','normal','fontsize',12)
 set(gca,'fontsize',12)
 grid on, box on
 set(gca,'fontsize',12)

 subplot(3,3,3)
 histogram(1e3*TRENDL,[0:.2:4],'normalization','probability','edgecolor','none','facecolor',clr0(2,:))
 grid on, box on, hold on
 histogram(1e3*TRENDL2,[0:.2:4],'normalization','probability','edgecolor','none','facecolor',clr0(4,:))
 xlabel('RSL Historical Trend (mm/yr)')
 ylabel('Probability')
 legend([{'CONUS'};{'No Western Gulf'}],'location','northwest','orientation','vertical'), legend boxoff
 axis([1 4 0 0.8])
 set(gca,'xtick',0:4,'ytick',0:0.2:0.8)
 set(gca,'fontsize',12)
 title([{'c. CONUS RSL'};{'Historical Trend (mm/yr)'}],'fontweight','normal')

 subplot(3,3,6)
 histogram(1e3*ACCQ,[0:0.002:0.04],'normalization','probability','edgecolor','none','facecolor',clr0(2,:))
 grid on, box on, hold on
 histogram(1e3*ACCQ2,[0:0.002:0.04],'normalization','probability','edgecolor','none','facecolor',clr0(4,:))
 xlabel('RSL Acceleration (mm/yr^2)')
 ylabel('Probability')
 legend([{'CONUS'};{'No Western Gulf'}],'location','northwest','orientation','vertical'), legend boxoff
 axis([0 0.04 0 0.4])
 set(gca,'xtick',0:0.01:0.04,'ytick',0:0.1:0.4)
 set(gca,'fontsize',12)
 title([{'d. CONUS RSL'};{'Acceleration (mm/yr^2)'}],'fontweight','normal')

 subplot(3,3,9)

 hold on, box on, grid on
 plot(1900:2024,1e3*prctile(RATEL,50),'color',clr0(3,:),'linewidth',2)
 c=fill([1900:2024 2024:-1:1900],1e3*[prctile(RATEL,2.5) fliplr(prctile(RATEL,97.5))],clr0(3,:));
 set(c,'edgecolor','none','FaceAlpha',0.2)
 plot(1900:2024,1e3*prctile(RATEQ,50),'color',clr0(2,:),'linewidth',2)
 plot(1900:2024,1e3*prctile(RATEQ,97.5),'--','color',clr0(2,:),'linewidth',1)
 plot(1900:2024,1e3*prctile(RATEQ,2.5),'--','color',clr0(2,:),'linewidth',1)
 legend([{'Historical Trend (Median)'};{'Historical Trend (95% CI)'};{'Quadratic Rate (Median)'};{'Quadratic Rate (95%CI)'}],'location','southeast','fontsize',7), legend boxoff
 axis([1900 2025 1 5])
ylabel('Rate (mm/yr)','fontsize',12,'fontweight','normal')
set(gca,'fontsize',12)
title('e. CONUS RSL Rate (mm/yr)','fontweight','normal')
set(gca,'fontsize',12,'fontweight','normal','xtick',1900:25:2025,'ytick',0:6)
set(gca,'fontsize',12)
xtickangle(45)
