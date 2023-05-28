function [SatECI,Satlla,DateVector] = FP01_CreateSatGeometry(startTime,stopTime,Param,Elem)

% The satelliteScenario object represents a 3D arena consisting of satellites, 
% ground stations, and the interactions between them. 
% Use this object to model satellite constellations, 
% model ground station networks, 
% perform access analyses between the satellites and the ground stations, 
% and visualize the results.
sc  = satelliteScenario(startTime,stopTime,Param.ts);

% adds a Satellite object from Keplerian elements defined in the Geocentric Celestial Reference Frame (GCRF) to the satellite scenario.
sat = satellite(sc,Elem.a,Elem.e,Elem.Inc,Elem.RAAN,Elem.omega,Elem.TA,'OrbitPropagator','two-body-keplerian','name','SAR sat');

% returns a 3-by-n-by-m array of the position history pos of each satellite in the vector sat, where n is the number of 
% time samples and m is the number of satellites. 
% The rows represent the x, y, and z coordinates of the satellite in the Geocentric Celestial Reference Frame (GCRF).
[SatECI,~,DateTime] = states(sat);

% converts the datetime or duration value t to a date vector—that is, 
% a numeric vector whose six elements represent the year, month, 
% day, hour, minute, and second components of t.
DateVector = datevec(DateTime);
DateVector(end,:)=[];       % Last sample is coming sometimes corrupted
SatECI(:,end)=[]; % Trim the last reading it has some errors. SatECI = Earth centered intertial

% converts Earth-centered inertial (ECI) coordinates, specified by position, to latitude, longitude, altitude (LLA) geodetic coordinates.
% The conversion is based on the Coordinated Universal Time (UTC) you specify.
Satlla = eci2lla(SatECI',DateVector);
%Satlla(:,3) = ones(size(Satlla,1),1)*Param.h; % Assume an ideal sphere, added 2023 to simplify the range migration process (it makes the approach distance to the center of the swath symetrical)
end