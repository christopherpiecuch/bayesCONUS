%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ClusMult = determine_clusters(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used in Piecuch 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is specific to the spatial covariance structure of sea level
% variation along the U.S. East Coast.  Given a set of N locations, this
% function creates an [NxN] matrix ("C") filled with zeros and ones.  The
% values C_ij equal 1 iff locations i and j are both *either* north or
% south of Cape Hatteras, but 0 if locations i and j are on opposite sites
% of hatteras.
% INPUT: 
%   Y           Latitudes of locations
% OUTPUT:
%   ClusMult    The "C" matrix defined above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by CGP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ClusMult = determine_clusters(X)

X=double(X>-110);
X(X==0)=-1;
ClusMult=X'*X;
ClusMult(ClusMult<0)=0;


return