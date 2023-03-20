function [Satlla,latSawthMid,lonSwathMid,slantrangeMid,Swathwidths_m,latSawthL1,lonSwathL1,slantrange1,latSawthL2,lonSwathL2,slantrange2,sataz]=F02_FindSwath(SatECI,DateVector,RadPar,E)
% This is to compute the approximate azimuth of the swath (also the
% satellite azimuth -> i.e. the direction of motion of the satellite

% Convert Earth-centered inertial (ECI) coordinates of satellite into latitude, longitude, altitude (LLA) geodetic coordinates
Satlla = eci2lla(SatECI',DateVector);    % The conversion is based on the Universal Coordinated Time (UTC) specified by Date vector

if RadPar.Left == 1
    % Adding 90deg if the scanning on the left side of the trajectory
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2)) +90;
else
    % Subtracting 90deg if the scanning on the rigth side of the trajectory
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2)) -90;
end

% Finding Mid-swath position (Longitude and latitude)
[latSawthMid,lonSwathMid,slantrangeMid] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),sataz,RadPar.AntOffNadir,E);

% Finding 2 edges of the swath (Longitude and latitude) along azimuth
[latSawthL1,lonSwathL1,slantrange1] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir-RadPar.SwathBeam/2,E);
[latSawthL2,lonSwathL2,slantrange2] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir+RadPar.SwathBeam/2,E);

% Finding the Swath width in meters
[Swathwidths_m,~] = distance(latSawthL1(1),lonSwathL1(1),latSawthL2(1),lonSwathL2(1));
Swathwidths_m = deg2km(Swathwidths_m)*1000;
end

