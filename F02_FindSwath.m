function [latSawthMid,lonSwathMid,slantrangeMid,Swathwidths_m,latSwathL1,lonSwathL1,latSwathL2,lonSwathL2,slantrange1,slantrange2,sataz]=F02_FindSwath(Satlla,RadPar,E)
% This is to compute the approximate azimuth of the swath (also the satellite azimuth -> i.e. the direction of motion of the satellite

% This is to compute the azimuth for each point of the satellite motion
for eta=1:size(Satlla,1)-10 
    if RadPar.Left == 0  % for the case from South to North - RadPar.Left == 1 for the case from North to South 
        % Adding 90deg if the scanning on the left side of the trajectory
        sataz(eta) = azimuth(Satlla(eta,1),Satlla(eta,2),Satlla(eta+1,1),Satlla(eta+1,2),E) +90;
    else
        % Subtracting 90deg if the scanning on the rigth side of the trajectory
        sataz(eta) = azimuth(Satlla(eta,1),Satlla(eta,2),Satlla(eta+1,1),Satlla(eta+1,2),E) -90;
    end
end

p = polyfit(1:size(Satlla,1)-10,sataz,1);
sataz = p(1) *(1:size(Satlla,1))+ p(2);
sataz=sataz';
 
% Finding Mid-swath position (Longitude and latitude)
[latSawthMid,lonSwathMid,slantrangeMid] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),sataz,RadPar.AntOffNadir,E);

% Finding 2 edges of the swath (Longitude and latitude) along azimuth
[latSwathL1,lonSwathL1,slantrange1] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir-RadPar.SwathWidthDeg/2,E);
[latSwathL2,lonSwathL2,slantrange2] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir+RadPar.SwathWidthDeg/2,E);

% Finding the Swath width in meters
[Swathwidths_m,~] = distance(latSwathL1(1),lonSwathL1(1),latSwathL2(1),lonSwathL2(1),E);
end