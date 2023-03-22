% Orbital Paramters / Elements
c=physconst('LightSpeed');              % Speed of light C [m/s]
ColorOrder =colororder;
E=wgs84Ellipsoid;                       % Reference Ellipsoide WGS84
%E=referenceSphere('Earth');            % Ideal sphere
h= 300e3;                               % Height of the platform [m]
Re = earthRadius;                       % Store Earth radius into Re variable
Elem.a = Re+h;                          % Semi-major axis of the platform
Elem.e   = 0;                           % Eccentricity - 0 means circular orbit, non zero for eliptical orbit 0 < e < 1
Elem.omega = 90;                        % Argument of periapsis in degrees - Defines the position of a body moving from vertical axis
Elem.Inc = 88;                          % Inclination in degrees - Angle between the orbital plane and the equator between 0 and 180 degrees
Elem.TA  = 123.78;                      % True Anomaly in degrees - Satalliete position in the orbit
Elem.RAAN = 192.52;                     % Right Ascension Of Ascending Node in degress - Ω angle is measured eastwards from x axis to ascending node on nodal line
Param.mu = 3.986e14;                    % Earth's standard gravitational parameter
Param.T = 2*pi*sqrt(Elem.a^3/Param.mu); % Orbital period
Param.dt = 0.003;                       % Orbit time step - Geo-sampling period in seconds - Slowtime sampling
Param.PRF = 1/Param.dt;                 % Pulse Repeatition Frequency [PRF]
Param.ft = 0.8;                         % flight duration - along azimuth direction in seconds
%% Targets
Param.NtargetsRange = 200;              % Number of targets in each eta bin
%% Time
startTime = datetime('01-Jan-2022 08:00:00');        % Set up start time
stopTime = startTime + seconds(Param.ft);            % Set up end time
%%
RadPar.fo = 300e6;                      % Carrier frequency [Hz] ==> Hypothetical parameters 
RadPar.Lambda = freq2wavelen(RadPar.fo);% Wavelength [m]
RadPar.AntOffNadir = 30;                % Antenna Offnadir angle [pointing angle] [degrees]
RadPar.BeamRange = 5;                   % Radar Beamwidth in the range direction [degrees]
RadPar.Ns = 1;                          % Number of sub-strips - N=0.5, BeamAz = 1.8466, N=1, BeamAz = 0.92, N=2, BeamAz = 0.46, N=3, BeamAz = 0.31
RadPar.BeamAz = 15;                     % Radar Beamwidth in the azimuth direction [degrees] - RadPar.BeamAz = 2*atand(Ls/(2*Ro)) , Ls =Param.ft*speed/RadPar.Ns 

RadPar.SwathWidthDeg = 0.5;             % Swath Beamwidth [degrees] - To define the swath lines of the targets
RadPar.Gain  = 20^(12/10);              % Antenna Gain 12 dBi example
RadPar.bw = 20e6;                       % Radar Bandwidth [Hz] ==> Hypothetical parameters
RadPar.fs = RadPar.bw*2;                % Sampling rate (Hz] - Fast time sampling frequency ( bw <~ fs < 10*bw) fs = alfa*Bw, alfa is between 1.1 and 2 (the multiplier is added to improve image quality)
RadPar.ts = (1/RadPar.fs);              % Sample time - 1/fs [s] - Fasttime sampling
RadPar.T = 2e-6;                        % Pulse width - Tp [s]
RadPar.K =  (RadPar.bw /RadPar.T);      % Chirp rate - K
RadPar.Left = 1;                        % The scanning is on the left side of the satellite trajectory
