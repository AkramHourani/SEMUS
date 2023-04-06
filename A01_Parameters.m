% Orbital Paramters / Elements
c=physconst('LightSpeed');                      % [m/s] Speed of light C 
ColorOrder =colororder;
E=referenceSphere('Earth');                     % Reference ideal sphere
Param.h= 600e3;                                 % [m] Height of the platform 
Re = earthRadius;                               % Store Earth radius into Re variable
Elem.a = Re+Param.h;                            % [m] Semi-major axis of the platform
Elem.e   = 0;                                   % Eccentricity - 0 means circular orbit, non zero for eliptical orbit 0 < e < 1
Elem.omega = 90;                                % [deg] Argument of periapsis - Defines the position of a body moving from vertical axis
Elem.Inc = 34.8;                                % [deg] Inclination - Angle between the orbital plane and the equator between 0 and 180 degrees
Elem.TA  = -175.72;                             % [deg] True Anomaly - Satalliete position in the orbit
Elem.RAAN = 90;                                 % [deg] Right Ascension Of Ascending Node - Ω angle is measured eastwards from x axis to ascending node on nodal line
Param.mu = 3.986e14;                            % [m3/s2] Earth's standard gravitational parameter
Param.T = 2*pi*sqrt(Elem.a^3/Param.mu);         % [s] Orbital period
Param.PRF = 2e3;                                % [Hz] Pulse Repeatition Frequency - PRF
Param.ts = 1/Param.PRF;                         % [s] Orbit time step - Geometric-sampling period - Slowtime sampling
Param.ScanDuration = seconds(1.5);              % [s] Flight duration - along azimuth direction
%% Targets
Param.NtargetsAz = 300;                         % Number of targets in each eta bin
Param.NtargetsRange = 300;                      % Number of targets in each range bin
Param.Margin = 1.2;                             % Range margin factor. This is a margin to include targets farther than swatch width, i.e. for larger squinet angle targets
%% Time
startTime = datetime('01-Jan-2022 08:00:00');   % [s] Set up the start time
stopTime  = startTime + Param.ScanDuration ;    % [s] Set up the end time
%% Radar Paramters
RadPar.fo = 1.3e9;                              % Carrier frequency [Hz] - L-band
% RadPar.fo = (500e6);                            % [Hz] Carrier frequency
RadPar.Lambda = freq2wavelen(RadPar.fo);        % [m] Wavelength
RadPar.AntOffNadir = 30;                        % [deg] Antenna Off-Nadir angle (pointing angle)

RadPar.BeamRange = 5;                           % [deg] Antenna beamwidth in the range direction
RadPar.BeamAz    = 5;                           % [deg] Antenna beamwidth in the azimuth direction
RadPar.Gain  = 10^(12/10);                      % [dBi] Antenna gain example

RadPar.bw = 15e6;                               % [Hz] Bandwidth of the RF signal to calcualte the ramp rate
RadPar.fs = RadPar.bw*2;                        % [Hz] Sampling rate (the multiplier is added to improve image quality)
RadPar.ts = (1/RadPar.fs);                      % [s] Sample time
RadPar.T = 5e-6;                                % [s] Pulse width
RadPar.K =  (RadPar.bw /RadPar.T);              % [Hz/s] Ramp (chirp) rate

RadPar.SwathWidthDeg = 0.5;                     % [deg] Swath width in degrees
RadPar.Left = 1;                                % The scanning on the left side of the satellite trajectory

