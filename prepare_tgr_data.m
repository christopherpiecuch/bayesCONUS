%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [DATA,X,Y,NAME] = prepare_tgr_data(la1,la2,lo1,lo2,minnum,coastcode,exceptions,time0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads annual sea level data from specified geographic region and
% coastline from the PSMSL database into into a [NxK] matrix, where N is
% the number of tide gauges and K is the number of years of the data
% record; lack of data are filled with NaNs.
% INPUT: 
%   la1         Southern latitudinal bounds of study region
%   la2         Northern latitudinal bounds "
%   lo1         Western longitudinal bounds "
%   lo2         Eastern longitudinal bounds "
%   minnum      Minimum number of data points to consider a gauge record
%   coastcode   PSMSL coastline ID  
%   exceptions  Additional tide gauge IDs to include (that might not meat
%               the coastcode or minnum requirements but include anyway)
%   time0       Initial time
% OUTPUT:
%   DATA        [NxK] array of sea level values in units of m
%   X           Longitude of N tide gauges
%   Y           Latitude of N tide gauges
%   NAME        Name of N tide gauge sites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DATA,X,Y,NAME,COAST] = prepare_tgr_data2(la1,la2,lo1,lo2,minnum,coastcode,time0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load annual tide gauge data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dircur=pwd;
cd rlr_annual_20250802/ % switch to PSMSL annual data directory
 % The above directory was downloaded from www.psmsl.org/data/ 
 % See PSMSL website for more recent data
 data = readAnnual([dircur,'/rlr_annual_20250802']);
cd(dircur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete empty entries in database or records not along specified
% coastlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete records outside the specified geographic region or with fewer than
% "minnum" data points.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=numel(data):-1:1
    %if ~ismember(data(n).id,exceptions)
%        if (sum(~isnan(data(n).height))<=minnum) 
        if ((max(data(n).year)-min(data(n).year))<=minnum) 
            data(n)=[];
        elseif (data(n).latitude<la1)||(data(n).latitude>la2)||(data(n).longitude<lo1)||(data(n).longitude>lo2)
            data(n)=[]; 
        end
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort remaining records by latitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:numel(data)
 lat2sort(n)=data(n).latitude;
end
[Y,I]=sort(lat2sort);
data2=data(I);
clear data
data=data2; clear data2 I Y lat2sort

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Be careful of flagged or interpolated values in PSMSL database. Make them
% into NaNs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:numel(data)
    iTakeOut=[]; iTakeOut=find((data(n).dataflag)|(data(n).interpolated));
    data(n).height(iTakeOut)=nan;
end

% remove data less than min num and not on coastline and doesn't have
% recent data
for kk=numel(data):-1:1
 testval=sum(~isnan(data(kk).height));
 latestyear=data(kk).year(find(~isnan(data(kk).height),1,'last'));
 if testval<minnum|~ismember(data(kk).coastline,coastcode)|latestyear<2020
  data(kk)=[];
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert data structure into matrix array structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=numel(data);
for n=1:N, %disp(num2str(n))
    if n==1
        T_0=data(n).year(1);
        T_F=data(n).year(numel(data(n).year));
    else
        if T_0>data(n).year(1)
           T_0=data(n).year(1);
        end
        if T_F<data(n).year(numel(data(n).year))
           T_F=data(n).year(numel(data(n).year));
        end
    end
end
K=T_F-T_0+1;
T=T_0:T_F;
DATA=nan(N,K);
for n=1:N
   t_0=find(T==data(n).year(1));
   t_f=find(T==data(n).year(numel(data(n).year)));
   DATA(n,t_0:t_f)=1e-3*data(n).height; % Convert data from mm to m
   X(n)=data(n).longitude;
   Y(n)=data(n).latitude;
   NAME(n).name=data(n).name;
   COAST(n)=data(n).coastline;
   clear t_0 t_f
end

i=find(T==1900);
if ~isempty(time0)
    i=find(T==time0);
end
DATA(:,1:(i-1))=[];
T(1:(i-1))=[];
T_0=T(1);
N=numel(data);

% adjust for any weird datums
for n=1:N
 if (nanmean(DATA(n,:))>7.2)||(nanmean(DATA(n,:))<6.8)
  DATA(n,:)=DATA(n,:)-nanmean(DATA(n,:))+7;
 end
end

return