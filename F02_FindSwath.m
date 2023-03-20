function [latSawthMid,lonSwathMid,slantrangeMid,Swathwidths_m,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,sataz]=F02_FindSwath(Satlla,RadPar,E)

% This is to compute the azimuth of the satellite motion
for eta=1:size(Satlla,1)-10 % Compute the azimuth for each point
if RadPar.Left == 1
    sataz(eta) = azimuth(Satlla(eta,1),Satlla(eta,2),Satlla(eta+1,1),Satlla(eta+1,2),E) +90;
else
    sataz(eta) = azimuth(Satlla(eta,1),Satlla(eta,2),Satlla(end+1,1),Satlla(eta+1,2),E) -90;
end
end

p = polyfit(1:size(Satlla,1)-10,sataz,1);
sataz = p(1) *(1:size(Satlla,1))+ p(2);
sataz=sataz';

[latSawthMid,lonSwathMid,slantrangeMid] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),sataz,RadPar.AntOffNadir,E);
[latSawthL1,lonSwathL1,~] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir-RadPar.SwathWidthDeg/2,E);
[latSawthL2,lonSwathL2,~] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir+RadPar.SwathWidthDeg/2,E);

%Find the Swath width in meters
[Swathwidths_m,~] = distance(latSawthL1(1),lonSwathL1(1),latSawthL2(1),lonSwathL2(1));
Swathwidths_m = deg2km(Swathwidths_m)*1000;
end