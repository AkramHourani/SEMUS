function [latSawthMid,lonSwathMid,slantrangeMid,Swathwidths_m,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,sataz]=F02_FindSwath(Satlla,RadPar,E)

% This is to compute the approximate azimuth of the swath (also the
% satellite azimuth -> i.e. the direction of motion of the satellite
if RadPar.Left == 1
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2)) +90;
else
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2)) -90;
end

[latSawthMid,lonSwathMid,slantrangeMid] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),sataz,RadPar.AntOffNadir,E);
[latSawthL1,lonSwathL1,~] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir-RadPar.SwathWidth/2,E);
[latSawthL2,lonSwathL2,~] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir+RadPar.SwathWidth/2,E);

%Find the Swath width in meters
[Swathwidths_m,~] = distance(latSawthL1(1),lonSwathL1(1),latSawthL2(1),lonSwathL2(1));
Swathwidths_m = deg2km(Swathwidths_m)*1000;
end