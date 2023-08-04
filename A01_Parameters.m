%% Orbital Paramters / Elements
c=physconst('LightSpeed');                      % [m/s] Speed of light C 
ColorOrder =colororder;
E=referenceSphere('Earth');                     % Reference ideal sphere
Param.h= 600e3;                                 % [m] Height of the platform 
Re = earthRadius;                               % Store Earth radius into Re variable
Elem.a = Re+Param.h;                            % [m] Semi-major axis of the platform
Elem.e   = 0;                                   % Eccentricity - 0 means circular orbit, non zero for eliptical orbit 0 < e < 1
Elem.omega = 90;                                % [deg] Argument of periapsis - Defines the position of a body moving from vertical axis
Elem.Inc = 97.44;                               % [deg] Inclination - Angle between the orbital plane and the equator between 0 and 180 degrees
Elem.TA  = -127.46;                             % [deg] True Anomaly - Satalliete position in the orbit
Elem.RAAN = -25.63;                             % [deg] Right Ascension Of Ascending Node - Ω angle is measured eastwards from x axis to ascending node on nodal line
Param.mu = 3.986e14;                            % [m3/s2] Earth's standard gravitational parameter
Param.Tf = 2*pi*sqrt(Elem.a^3/Param.mu);        % [s] Orbital period
Param.PRF = 2e3;                                % [Hz] Pulse Repeatition Frequency - PRF
Param.tg = 1/Param.PRF;                         % [s] Orbit time step - Geometric-sampling period - Slowtime sampling
%% Timing Parameters
Param.ScanDuration = seconds(1);                % [s] Flight duration - along azimuth direction
startTime = datetime('01-Jan-2023 6:00:00');    % [s] Set up the start time
stopTime  = startTime + Param.ScanDuration ;    % [s] Set up the end time
%% Targets' Parameters
Param.NtargetsAz = 500;                         % Number of targets in each eta bin
Param.NtargetsRange = 500;                      % Number of targets in each range bin
Param.Margin = 1.2;                             % Range margin factor. This is a margin to include targets farther than swatch width, i.e. for larger squinet angle targets
%% Radar and Antenna Paramters
RadPar.fo = 1.2e9;                              % [Hz] Carrier frequency - L-band
RadPar.Lambda = freq2wavelen(RadPar.fo);        % [m] Wavelength
RadPar.bw = 20e6;                               % [Hz] Bandwidth of the RF signal to calcualte the ramp rate
RadPar.fs = RadPar.bw*2;                        % [Hz] Sampling rate (the multiplier is added to improve image quality)
RadPar.ts = (1/RadPar.fs);                      % [s] Sample time
RadPar.T = 5e-6;                                % [s] Pulse width
RadPar.K = (RadPar.bw /RadPar.T);               % [Hz/s] Ramp (chirp) rate
RadPar.AntOffNadir = 35;                        % [deg] Antenna Off-Nadir angle (pointing angle)
RadPar.BeamRange = 4;                           % [deg] Antenna beamwidth in the range direction
RadPar.BeamAz  = 0.2;                             % [deg] Antenna beamwidth in the azimuth direction
RadPar.SwathWidthDeg = 0.5;                     % [deg] Swath width in degrees
RadPar.Left = 1;                                % The scanning on the left side of the satellite trajectory
%% Power Parameters
RadPar.Gain  = 10^(12/10);                      % [dBi] Antenna max gain example
RadPar.Pt = 10^(35/10);                         % [dBW] Transmitted power example
RadPar.NESO = -24.5;                            % Noise equivalent sigma naught [dB] 
RadPar.FeederL = 10;                            % Feeder loss [dB]
RadPar.Duty = 0.1;                              % Duty Cycle
RadPar.kb = 1.38e-23;                           % Boltzmann's constant [J K−1]
RadPar.AT = 290;                                 % Absolute temperature [degrees Kelvin]
RadPar.Bn = Param.PRF/RadPar.Duty;              % Noise equivalent bandwidth [Hz]
RadPar.F = 8 ;                                  % Noise figure [dB]