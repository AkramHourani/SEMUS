function [SatECI,velocity,DateTime] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem)

% Create the satellite scenario and the orbit
sc = satelliteScenario(startTime,stopTime,Param.ts);
sat=satellite(sc,Elem.a,Elem.e,Elem.Inc,Elem.RAAN,Elem.omega,Elem.TA,'OrbitPropagator','two-body-keplerian',...
    'name','SAR sat');

% Extract the Geocentric Celestial Reference Frame (GCRF) | Earth-centered inertial (ECI) - Satellite coordinates and Datetime in UTC  
[SatECI,velocity,DateTime] = states(sat);
end