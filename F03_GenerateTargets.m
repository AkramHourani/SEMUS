function [Targetlat,Targetlon]= F03_GenerateTargets(latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,Param)
% This loop will walk step by step along the track (azimuth) to generate the sampling
% points (arc) between the two edges of the swath based on the great cricle
% connecting the two edges.

Targetlat=zeros(length(latSawthL1),Param.NtargetsRange+1);
Targetlon=zeros(length(latSawthL1),Param.NtargetsRange+1);
etaTotal=length(latSawthL1);

for eta=1:etaTotal
[Targetlat(eta,:),Targetlon(eta,:)] = gcwaypts(latSawthL1(eta),lonSwathL1(eta),latSawthL2(eta),lonSwathL2(eta),Param.NtargetsRange); 
end