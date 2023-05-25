function [Targetlat,Targetlon]= F04_GenerateTargets(latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,Param)
% This loop will walk step by step along the track (azimuth) to generate the sampling
% points (arc) between the two edges of the swath based on the great cricle
% connecting the two edges.

% Creat the variable tamplate
Targetlat=zeros(Param.NtargetsAz,Param.NtargetsRange);
Targetlon=zeros(Param.NtargetsAz,Param.NtargetsRange);


%Interpolate the swath lines to match the number of azimuth samples
[TargetEdgeLat1, TargetEdgeLon1] = gcwaypts(latSawthL1(1),lonSwathL1(1),latSawthL1(end),lonSwathL1(end),Param.NtargetsAz-1); 
[TargetEdgeLat2, TargetEdgeLon2] = gcwaypts(latSawthL2(1),lonSwathL2(1),latSawthL2(end),lonSwathL2(end),Param.NtargetsAz-1); 

for AzTarget=1:Param.NtargetsAz

[Targetlat(AzTarget,:),Targetlon(AzTarget,:)] = gcwaypts(TargetEdgeLat1(AzTarget),TargetEdgeLon1(AzTarget),TargetEdgeLat2(AzTarget),TargetEdgeLon2(AzTarget),Param.NtargetsRange-1); 
end