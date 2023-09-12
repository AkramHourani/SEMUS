function [latSawthMid,lonSwathMid,slantrangeMid,Swathwidths_m,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,sataz]=FP02_FindSwath(Satlla,RadPar,E)
% This is to compute the approximate azimuth of the swath (also the satellite azimuth -> i.e. the direction of motion of the satellite

% This is to compute the azimuth for each point of the satellite motion
for eta=1:size(Satlla,1)-10 
    if RadPar.Left == 1
        % Adding 90deg if the scanning on the left side of the trajectory
        sataz(eta) = azimuth(Satlla(eta,1),Satlla(eta,2),Satlla(eta+1,1),Satlla(eta+1,2),E) +90;
    else
        % Subtracting 90deg if the scanning on the rigth side of the trajectory
        sataz(eta) = azimuth(Satlla(eta,1),Satlla(eta,2),Satlla(eta+1,1),Satlla(eta+1,2),E) -90;
    end
end

% p = polyfit(x,y,n) returns the coefficients for a polynomial p(x) of degree n that is a best fit (in a least-squares sense) 
% for the data in y. The coefficients in p are in descending powers, and the length of p is n+1
p = polyfit(1:size(Satlla,1)-10,sataz,1);
sataz = p(1) *(1:size(Satlla,1))+ p(2);
sataz = sataz';
 
% [lat,lon,slantrange] = lookAtSpheroid(lat0,lon0,h0,az,tilt,spheroid) computes the latitude and longitude (lat and lon) 
% of the first point of intersection of the line-of-sight from a viewpoint in space with the surface of an oblate spheroid. 
% If the line-of-sight does not intersect with the spheroid, lat and lon contain NaNs.

% Finding Mid-swath position (Longitude and latitude)
[latSawthMid,lonSwathMid,slantrangeMid] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),sataz,RadPar.AntOffNadir,E);

% Finding 2 edges of the swath (Longitude and latitude) along azimuth
[latSawthL1,lonSwathL1,~] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir-RadPar.SwathWidthDeg/2,E);
[latSawthL2,lonSwathL2,~] = lookAtSpheroid(Satlla(:,1),Satlla(:,2),Satlla(:,3),...
    sataz,RadPar.AntOffNadir+RadPar.SwathWidthDeg/2,E);

% Finding the Swath width in meters
[Swathwidths_m,~] = distance(latSawthL1(1),lonSwathL1(1),latSawthL2(1),lonSwathL2(1),E);
end