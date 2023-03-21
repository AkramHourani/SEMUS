function [R,Ro,GRP,SAR_Dist_Edge1,SAR_Dist_Edge2,Swathwidth_SARDistance] = F03_FindGRP(latSawthMid,lonSwathMid,Satlla,E,Idx)

[~,~,R] = geodetic2aer(latSawthMid(Idx),lonSwathMid(Idx),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
% Find Reference range at the centre of the dwell at ground refernece point(GRP) - Find Index of mid point of the dwell
[Ro,~] = min(R);
% Swath center point (Ground Reference Point GRP)
GRP = [latSawthMid(Idx),lonSwathMid(Idx),0];
[~,~,SAR_Dist_Edge1] = geodetic2aer(latSawthL1(Idx),lonSwathL1(Idx),0,Satlla(Idx,1),Satlla(Idx,2),Satlla(Idx,3),E);
[~,~,SAR_Dist_Edge2] = geodetic2aer(latSawthL2(Idx),lonSwathL2(Idx),0,Satlla(Idx,1),Satlla(Idx,2),Satlla(Idx,3),E);
Swathwidth_SARDistance = abs(SAR_Dist_Edge1-SAR_Dist_Edge2);
end