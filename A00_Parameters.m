% Orbital Paramters
c=physconst('LightSpeed');
ColorOrder =colororder;
E=wgs84Ellipsoid;
%E=referenceSphere('Earth'); % Ideal sphere
h= 300e3;
Re = earthRadius;
Elem.a = Re+h; % Semi-major axis
Elem.e   = 0; % Eccentricity
Elem.omega = 90; % [Degrees] Argument of periapsis
Elem.omega = 90; % [Degrees] Argument of periapsis
Elem.Inc = 88;  % Inclination in degrees
Elem.TA  = 123.78;  %  in degrees
Elem.RAAN = 192.52;
Param.mu = 3.986e14; % Erth's standard gravitational parameter
Param.T = 2*pi*sqrt(Elem.a^3/Param.mu); % Orbital period
Param.dt = 0.003; % Orbit time step
%% Targets
Param.NtargetsRange = 200; % number of targets in each eta bin
%% Time
startTime = datetime('01-Jan-2022 08:00:00');
stopTime = startTime + seconds(0.8);
%%
RadPar.fo = (300e6); %Carrier frequency
RadPar.AntOffNadir = 30; % Antenna Offnadir angle (pointing angle)

RadPar.BeamRange = 5; %Beamwidth in the range diration
RadPar.BeamAz = 15; %Beamwidth in the azimuth diration
RadPar.SwathWidthDeg = 0.5; % Swath width in degrees
RadPar.Gain  = 20^(12/10); % antenna gain 12 dBi example
RadPar.bw = 20e6; % Bandwidth of the RF signal to calcualte the ramp rate
RadPar.fs = RadPar.bw*2; % Sampling rate (the multiplier is added to improve image quality)
RadPar.ts = (1/RadPar.fs); % Sample time
RadPar.T = 2e-6; % Pulse width
RadPar.K =  (RadPar.bw /RadPar.T); % ramp rate
RadPar.Left = 1; % the scanning on the left side of the satellite trajectory

