function [SatECI] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem)
sc = satelliteScenario(startTime,stopTime,Param.dt);
sat=satellite(sc,Elem.a,Elem.e,Elem.Inc,Elem.RAAN,Elem.omega,Elem.TA,'OrbitPropagator','two-body-keplerian',...
    'name','SAR sat');
SatECI = states(sat); % Extract the GCRF