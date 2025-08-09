%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [G_LON,G_LAT] = grab_coastal_cells(la1,la2,lo1,lo2,t_lon,...
%   t_lat,g_lon,g_lat,gridRes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code populates a set of longitudes and latitudes that constitute a
% regular target grid of coastal points
% INPUT: 
%   la1         Southern latitudinal bounds of study region
%   la2         Northern latitudinal bounds "
%   lo1         Western longitudinal bounds "
%   lo2         Eastern longitudinal bounds "
%   t_lon       Longitudes of tide gauges
%   t_lat       Latitudes of tide gauges
%   g_lon       Longitudes of GPS stations
%   g_lat       Latitudes of GPS stations
%   gridRes     Resolution of target grid
% OUTPUT:
%   G_LON       Vector of regular longitude points
%   G_LAT       Vector of regular latitude points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G_LON,G_LAT] = grab_coastal_cells(la1,la2,lo1,lo2,t_lon,...
    t_lat,g_lon,g_lat,gridRes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRID_LON = -180:gridRes:180;
GRID_LAT = -90:gridRes:90;
[GRID_LAT GRID_LON]=meshgrid(GRID_LAT,GRID_LON);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load coast data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load coastlines.mat
lat=coastlat; long=coastlon; clear coastl*
I=[]; I=find(lat<la1|lat>la2|long<lo1|long>lo2);
long(I)=[]; lat(I)=[];

%Destim = EarthDistances([[long' t_lon g_lon]',[lat' t_lat g_lat]']);
%Destim = Destim(1:numel(long),(numel(long)+1):(numel(long)+numel(t_lon)+numel(g_lon)));
%Destim = min(Destim');
%long(find(Destim>110))=[];
%lat(find(Destim>110))=[];

% now see which lat/lon boxes have coastal points
GRID_LAT_S=GRID_LAT-gridRes/2;
GRID_LAT_N=GRID_LAT+gridRes/2;
GRID_LON_W=GRID_LON-gridRes/2;
GRID_LON_E=GRID_LON+gridRes/2;

hasCoastPoint=size(GRID_LAT);
for nx=1:(360/gridRes+1), %disp(num2str(nx))
 for ny=1:(180/gridRes+1)
  xv=[]; xv=[GRID_LON_W(nx,ny) GRID_LON_E(nx,ny) GRID_LON_E(nx,ny) GRID_LON_W(nx,ny) GRID_LON_W(nx,ny)];
  yv=[]; yv=[GRID_LAT_S(nx,ny) GRID_LAT_S(nx,ny) GRID_LAT_N(nx,ny) GRID_LAT_N(nx,ny) GRID_LAT_S(nx,ny)];
  IN=[]; IN = inpolygon(long,lat,xv,yv);
  hasCoastPoint(nx,ny)=sum(IN);
 end
end
hasCoastPoint=double(hasCoastPoint~=0);

targetLat=GRID_LAT(find(hasCoastPoint));
targetLon=GRID_LON(find(hasCoastPoint));

% define output
G_LAT = targetLat;
G_LON = targetLon;
ijk=find(G_LON<=-90&G_LAT<=17);
G_LAT(ijk)=[];
G_LON(ijk)=[];

% now omit non-contiguous us coastlines
states = shaperead('usastatelo', 'UseGeoCoords', true);
states(11)=[]; % delete hawaii
states(2)=[]; % delete alaska
% now create lat lons
usalon=[]; usalat=[];
for kk=1:numel(states);
    usalon=[usalon states(kk).Lon];
    usalat=[usalat states(kk).Lat];
end
for kk=numel(G_LON):-1:1, disp(num2str(kk))
 dd=[]; dd=distance(G_LAT(kk),G_LON(kk),usalat,usalon,6371);
 if min(dd)>50
  G_LON(kk)=[]; G_LAT(kk)=[];
 end
end
% now delete great lakes & st lawrence
dd=[]; dd=find(G_LON>-100&G_LON<-75&G_LAT>=41);
G_LAT(dd)=[]; G_LON(dd)=[];
dd=[]; dd=find(G_LON>-100&G_LAT>=46);
G_LAT(dd)=[]; G_LON(dd)=[];







return
