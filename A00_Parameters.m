% Orbital Paramters
c=physconst('LightSpeed');
ColorOrder =colororder;
E=referenceSphere('Earth');     % Ideal sphere
Param.h= 300e3;
Re = earthRadius;
Elem.a = Re+Param.h;                  % Semi-major axis
Elem.e   = 0; % Eccentricity
Elem.omega = 90; % [Degrees] Argument of periapsis
Elem.Inc = 88;  % Inclination in degrees
Elem.TA  = 123.78;  %  in degrees
Elem.RAAN = 192.52;
Param.mu = 3.986e14; % Erth's standard gravitational parameter
Param.T = 2*pi*sqrt(Elem.a^3/Param.mu); % Orbital period
Param.PRF = 2e3;  % [Hz] Pulse repition freqnecy 
Param.ts = 1/Param.PRF; % Slow time, time-step also the goemtric sampling time step
Param.ScanDuration = seconds(0.8);
%% Targets
Param.NtargetsAz = 20; % number of targets in each eta bin
Param.NtargetsRange = 20; % number of targets in each eta bin
Param.Margin = 1.2; % Range margin factor. This is a margin to include targets farther than swatch width, i.e. for larger squinet angle targets

%% Time
startTime = datetime('01-Jan-2022 08:00:00');
stopTime  = startTime + Param.ScanDuration ;
%% Radar Paramters
RadPar.fo = (1000e6); %Carrier frequency
RadPar.AntOffNadir = 30; % Antenna Offnadir angle (pointing angle)

RadPar.BeamRange = 5;      % [deg] Antenna beamwidth in the range diration
RadPar.BeamAz    = 5;     % [deg] Antenna beamwidth in the azimuth diration
RadPar.Gain  = 10^(12/10); % [dBi ]Antenna gain example

RadPar.bw = 15e6;       % Bandwidth of the RF signal to calcualte the ramp rate
RadPar.fs = RadPar.bw*2; % Sampling rate (the multiplier is added to improve image quality)
RadPar.ts = (1/RadPar.fs); % Sample time
RadPar.T = 5e-6; % Pulse width
RadPar.K =  (RadPar.bw /RadPar.T); % ramp rate

RadPar.SwathWidthDeg = 0.5; % Swath width in degrees
RadPar.Left = 1; % the scanning on the left side of the satellite trajectory

