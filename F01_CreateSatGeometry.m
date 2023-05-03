function [SatECI,Satlla,DateVector] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem)
sc  = satelliteScenario(startTime,stopTime,Param.ts);
sat = satellite(sc,Elem.a,Elem.e,Elem.Inc,Elem.RAAN,Elem.omega,Elem.TA,'OrbitPropagator','two-body-keplerian',...
    'name','SAR sat');
[SatECI,~,DateTime] = states(sat);
DateVector = datevec(DateTime);
DateVector(end,:)=[];       % Last sample is comming sometimes corrupted
SatECI(:,end)=[];           % Trim the last reading it has some errors
Satlla = eci2lla(SatECI',DateVector);
Satlla(:,3) = ones(size(Satlla,1),1)*Param.h; % Assume an ideal sphere, added 2023 to simplify the range migration process (it makes the approach distance to the center of the swath symetrical)
end