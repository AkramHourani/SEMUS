% Orbital Paramters / Elements
clc; clear; close all hidden
c = physconst('LightSpeed');            % Speed of light C [m/s]
ColorOrder = colororder;
E = wgs84Ellipsoid;                     % Reference Ellipsoide WGS84
h = 600e3;                              % Height of the platform [m]
Re = earthRadius;                       % Store Earth radius into Re variable
Elem.a = Re+h;                          % Semi-major axis of the platform
Elem.e   = 0;                           % Eccentricity - 0 means circular orbit, non zero for eliptical orbit 0 < e < 1
Elem.omega = 90;                        % Argument of periapsis in degrees - Defines the position of a body moving from vertical axis
Elem.Inc = 34.8;                        % Inclination in degrees - Angle between the orbital plane and the equator between 0 and 180 degrees
Elem.TA  = -175.72;                     % True Anomaly in degrees - Satalliete position in the orbit
Elem.RAAN = 90;                         % Right Ascension Of Ascending Node in degress - Ω angle is measured eastwards from x axis to ascending node on nodal line
Param.mu = 3.986e14;                    % Earth's standard gravitational parameter
Param.T = 2*pi*sqrt(Elem.a^3/Param.mu); % Orbital period
Param.PRF = 200;                        % Pulse Repeatition Frequency [PRF]
Param.ts = 1/Param.PRF;                 % Orbit time step - Geo-sampling period in seconds - Slowtime sampling
Param.ft = 1.5;                         % flight duration - along azimuth direction in seconds
%% Targets
Param.NtargetsRange = 400;              % Number of targets in each eta bin
%% Time
startTime = datetime('01-Jan-2022 08:00:00');        % Set up start time
stopTime = startTime + seconds(Param.ft);            % Set up end time
%% Radar Parameters
% RadPar.fo = 1.3e9;                    % Carrier frequency [Hz] - L-band
% RadPar.fo = 5.9e9;                    % Carrier frequency [Hz] - C-band
RadPar.fo = 300e6;                      % Carrier frequency [Hz] ==> Hypothetical parameters 
RadPar.Lambda = freq2wavelen(RadPar.fo);% Wavelength [m]
RadPar.AntOffNadir = 30;                % Antenna Offnadir angle [pointing angle] [degrees]
RadPar.BeamRange = 5;                   % Radar Beamwidth in the range direction [degrees]
RadPar.Ns = 1;                          % Number of sub-strips - N=0.5, BeamAz = 1.8466, N=1, BeamAz = 0.92, N=2, BeamAz = 0.46, N=3, BeamAz = 0.31
RadPar.BeamAz = 0.92;                   % Radar Beamwidth in the azimuth direction [degrees] - RadPar.BeamAz = 2*atand(Ls/(2*Ro)) , Ls =Param.ft*speed/RadPar.Ns 

RadPar.SwathBeam = 0.5;                 % Swath Beamwidth [degrees] - To define the swath lines of the targets
RadPar.Gain  = 10^(12/10);              % Antenna Gain 12 dBi example
RadPar.bw = 10e6;                       % Radar Bandwidth [Hz] ==> Hypothetical parameters
RadPar.fs = RadPar.bw*1.2;              % Sampling rate (Hz] - Fast time sampling frequency ( bw <~ fs < 10*bw) fs = alfa*Bw, alfa is between 1.1 and 1.2
RadPar.ts = 1/RadPar.fs;                % Sample time - 1/fs [s] - Fasttime sampling
RadPar.T = 10e-6;                       % Pulse width - Tp [s]
RadPar.K = RadPar.bw /RadPar.T;         % Chirp rate - K
RadPar.Left = 1;                        % The scanning on the left side of the satellite trajectory
%% Noise Parameters
Noise.Pt = 2400;                        % Assumed Transmitted power [Watt]
Noise.NESO = -24.5;                     % Noise equivalent sigma naught [dB] 
Noise.FeederL = 10;                     % Feeder loss [dB]
Noise.Duty = 0.1;                       % Duty Cycle
Noise.kb = 1.38e-23;                    % Boltzmann's constant [ J K−1]
Noise.T = 290;                          % Absolute temperature in degrees Kelvin
Noise.Bn = Param.PRF/Noise.Duty;        % Noise equivalent bandwidth
Noise.F = 8 ;                           % Noise figure [dB]
Noise.SNR = 24;                         % Assumed Signal to noise ratio  [dB]
%% AM Signal Parameters
AM.fc = RadPar.fo;                      % Carrier frequency same as SAR [Hz]
AM.fs = 50*RadPar.fs;                   % Sampling frequency same higher than SAR sampling frequency [Hz]
AM.t = (0:1/AM.fs:0.2);                 % Time base vector for the carrier modulated signal
AM.fm = 1e6;                            % Baseband signal frequency [Hz]       
AM.NumberofAM = 1;                      % Number of AM signals
AM.DeltaTheta = 0.05;                      % Shift of AM Transmitter from GRP longitude
AM.Gain = 1;                            % Asumme omin-directional isotropic AM transmitter G = 1
AM.SIR = 20;                            % Assumed SIR in [dB]
%% QPSK Signal Parameters
QPSK.NofBits=2^14;                      % Number of Transmitted bits
QPSK.fc = RadPar.fo;                    % Carrier frequency same as SAR sampling frequency [Hz]
QPSK.OF = 100;                          % Oversampling factor (multiples of fc) - at least 4 is better
QPSK.fs = RadPar.fs;                    % Sampling frequency same as SAR [Hz]
QPSK.NumberofQPSK = 1;                  % Number of QPSK signals
QPSK.DeltaTheta = 0.05;                    % Shift of QPSK Transmitter from GRP longitude
QPSK.Gain = 1;                          % Asumme omin-directional isotropic QPSK transmitter G = 1
QPSK.SIR = 22;                          % Assumed SIR in [dB]
%% LoRa Signal Parameters
% Generate 0.6595s of LORA data ==> Controlled by the SF=11
% SF ∈ {7, 8, 9, 10, 11, 12}. % Generate 0.16s of LORA data ==> Controlled by the SF=9
LORA.SF = 9;                            % Spreading factor
LORA.BW = 125e3 ;                       % Signal bandwidth of LoRa transmission [Hz]
% LORA.fc = RadPar.fo + 1e6 ;           % Carrier frequency [Hz]
LORA.fs = RadPar.fs ;                   % Sampling frequency-Same as radar range sampling frequency [Hz]
LORA.Delta_f = 1e6 ;                    % Frequency Shift [Hz]
LORA.NumberofLoRa = 1;                  % Number of LoRa signals
LORA.DeltaTheta = 0.05;                    % Shift of LoRa Transmitter from GRP longitude
LORA.Gain = 1;                          % Asumme omin-directional isotropic LORA transmitter G = 1
LORA.SIR = 22;                          % Assumed SIR in [dB]




