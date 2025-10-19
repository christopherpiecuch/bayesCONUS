%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function bayes_plot_residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 19 October 2025
% This code plots residual differences between
% Bayesian model solutions and raw tide-gauge data
% at all tide-gauge locations. The residuals should 
% look like white noise with comparable variance
% across all site (model assumption is that residuals
% are iid white noise in space and time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bayes_plot_residuals(experimentName,runnum)

% load files
a=[]; b=[]; bg=[]; l=[]; y=[]; yg=[];
for kk=1:numel(runnum), disp(num2str(kk))
 load(['bayes_model_solutions/experiment_',experimentName,'_runNum_',num2str(runnum(kk)),'.mat']);
 a=[a; A];
 l=[l; ELL];
 y=[y; Y];
 clear A ELL Y YG
 DATA=TGR_DATA;
end
A=a; L=l; Y=y; clear a b bg l y yg
NumIter=max(runnum)*numel(R);
T=1900:2024; T=T-mean(T);

% predict data from model solution
LD=permute(reshape(reshape(L,NumIter*numel(NAME),1)*ones(1,numel(T)),NumIter,70,125),[1 3 2]);
AD=permute(reshape(reshape(A,NumIter*numel(NAME),1)*T,NumIter,70,125),[1 3 2]);
YD=Y+LD+AD;

% plot tide gauge residuals (should all be white noise with comparable standard deviations)
sites=1:numel(LON);
sites=reshape(sites,5,numel(LON)/5);
[XS YS]=size(sites);
for kk=1:YS
 fig=figure('color','white');
 fig.Position(4) = fig.Position(4)*2;
 for ll=1:5
  subplot(5,1,ll)
  plot(1900:2024,prctile(squeeze(YD(:,:,sites(ll,kk))),50)-DATA(sites(ll,kk),:),'k','linewidth',2)
  hold on, box on, grid on
  plot(1900:2024,prctile(squeeze(YD(:,:,sites(ll,kk))),2.5)-DATA(sites(ll,kk),:),'k--','linewidth',1)
  plot(1900:2024,prctile(squeeze(YD(:,:,sites(ll,kk))),97.5)-DATA(sites(ll,kk),:),'k--','linewidth',1)
  title(NAME(sites(ll,kk)).name,'fontweight','normal')
  axis([1900 2025 -.06 .06])
  set(gca,'ytick',-.06:.03:.06)
  ylabel('Residuals (m)')
  legend([{'Median'};{'95% CI'}],'location','northwest','orientation','horizontal'), legend boxoff
 end
end