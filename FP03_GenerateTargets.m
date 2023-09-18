function [Targetlat,Targetlon]= FP03_GenerateTargets(latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,Param)

% This loop will walk step by step along the track (azimuth) to generate the sampling
% points (arc) between the two edges of the swath based on the great cricle
% connecting the two edges.
% Called once in ASP01_GenRawBi, the one that generates target for the
% Passive SAR
% The input is ARRAY of latSawthL1,lonSwathL1,latSawthL2,lonSwathL2

% Create the variable template
% For example : 
% Param.NtargetsAz = 11
% Param.NtargetsRange = 11

Targetlat = zeros(Param.NtargetsAz,Param.NtargetsRange);
Targetlon = zeros(Param.NtargetsAz,Param.NtargetsRange);

% Function : [lat,lon] = gcwaypts(lat1,lon1,lat2,lon2) returns the coordinates of equally 
% spaced points along a great circle path connecting two endpoints, (lat1,lon1) and (lat2,lon2).

% Interpolate the swath lines to match the number of azimuth samples
% For example for Active SAR (Illuminator of opportunity) : 
%       ----------------------------------------- L1
%
%
%       ----------------------------------------- L2
%
% For L1 :
% latSwathL1SoI(1) =  -33.7873
% lonSwathL1SoI(1) = 151.3131
% latSawthL1SoI(end) = -33.8381
% lonSwathL1SoI(end) = 151.2428 
%
% For L2 : 
% latSwathL2SoI(1) =  -33.7848
% lonSwathL2SoI(1) = 151.2457
% latSawthL2SoI(end) = -33.8356
% lonSwathL2I(end) = 151.2426 

% The results :
% TargetEdgeLat1 =      TargetEdgeLon1 =
%   -33.7873            151.3131
%   -33.7924            151.3128
%   -33.7975            151.3126
%   -33.8026            151.3123
%   -33.8076            151.3120
%   -33.8127            151.3117
%   -33.8178            151.3114
%   -33.8229            151.3111
%   -33.8279            151.3108
%   -33.8330            151.3106
%   -33.8381            151.3103
% 
% TargetEdgeLat2 =      TargetEdgeLon2 =
%   -33.7848            151.2457
%   -33.7899            151.2454
%   -33.7950            151.2452
%   -33.8000            151.2449
%   -33.8051            151.2446
%   -33.8102            151.2443
%   -33.8153            151.2440
%   -33.8203            151.2437
%   -33.8254            151.2434
%   -33.8305            151.2431
%   -33.8356            151.2428
%   
% Now we have : 
%      This is edge -->  o---o---o---o---o---o---o---o---o---o---o L1 
%
%
%      This is edge -->  o---o---o---o---o---o---o---o---o---o---o L2
%
% Similarly, calculate gcwaypts between point (pairs) 
%      This is edge -->   o---o---o---o---o---o---o---o---o---o---o L1 
%                         |   |   |   |   |   |   |   |   |   |   |
%                         o   o   o   o   o   o   o   o   o   o   o
%                         .   .   .   .   .   .   .   .   .   .   .
%                         o   o   o   o   o   X   o   o   o   o   o
%                         .   .   .   .   .   .   .   .   .   .   .
%                         o   o   o   o   o   o   o   o   o   o   o
%                         |   |   |   |   |   |   |   |   |   |   |
%     This is edge -->    o---o---o---o---o---o---o---o---o---o---o L2
%
%   TargetLat 11 x 11 ; TargetLon 11 x 11
% if GRP is chosen then the result will be [Targetlat,Targetlon] = [-33.8114,151.2781]
%
%
%
[TargetEdgeLat1, TargetEdgeLon1] = gcwaypts(latSawthL1(1),lonSwathL1(1),latSawthL1(end),lonSwathL1(end),Param.NtargetsAz-1); 
[TargetEdgeLat2, TargetEdgeLon2] = gcwaypts(latSawthL2(1),lonSwathL2(1),latSawthL2(end),lonSwathL2(end),Param.NtargetsAz-1); 

for AzTarget=1:Param.NtargetsAz
    [Targetlat(AzTarget,:),Targetlon(AzTarget,:)] = gcwaypts(TargetEdgeLat1(AzTarget),TargetEdgeLon1(AzTarget),TargetEdgeLat2(AzTarget),TargetEdgeLon2(AzTarget),Param.NtargetsRange-1); 
end

end