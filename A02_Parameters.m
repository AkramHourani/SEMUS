%% Noise and Interferers Parameters
%% AWGN signal parameter
Noise.SNR = 60;                         % Assumed Signal to noise ratio  [dB]
%% LoRa Signal Parameters
% Generate 0.6595s of LORA data ==> Controlled by the SF=11
% SF ∈ {7, 8, 9, 10, 11, 12}.           % Generate 0.16s of LORA data ==> Controlled by the SF=9
LORA.SF = 9;                           % Spreading factor
LORA.BW = 0.08e6 ;                       % Signal bandwidth of LoRa transmission [Hz]
% LORA.fc = RadPar.fo + 1e6 ;           % Carrier frequency [Hz]
LORA.Delta_f = 1e6 ;                    % Frequency Shift [Hz]
LORA.NumberofLoRa = 1;                  % Number of LoRa signals
LORA.latShift = 0.01;                   % Shift of LoRa Transmitter from GRP latitude
LORA.lonShift = -0.028;                 % Shift of LoRa Transmitter from GRP longitude
LORA.Gain = 1;                          % Asumme omin-directional isotropic LORA transmitter G = 1
LORA.SIR = -2;                          % Assumed SIR in [dB]
%% AM Signal Parameters
AM.fc = 1.01*RadPar.fo;                      % Carrier frequency same as SAR [Hz]
AM.BW = 0.01e6 ;                         % Signal bandwidth of LoRa transmission [Hz]
AM.t  = 0: 1/RadPar.fs :time2num(Param.ScanDuration); % Time base vector for the carrier modulated signal
AM.fm = 1.2*RadPar.fo;                           % Modulation frequency [Hz]       
AM.Index = 0.8;                      % Number of AM signals
AM.latShift = 0.01;                     % Shift of AM Transmitter from GRP latitude
AM.lonShift = -0.028;                     % Shift of AM Transmitter from GRP longitude
AM.Gain = 1;                            % Asumme omin-directional isotropic AM transmitter G = 1
AM.SIR = -20;                            % Assumed SIR in [dB]
%% QPSK Signal Parameters
QPSK.SymbolRate = 1e3;                  % Symbol rate
QPSK.fc = 1.01*RadPar.fo;               % Carrier frequency same as SAR sampling frequency [Hz]
QPSK.OF = 100;                          % Oversampling factor (multiples of fc) - at least 4 is better
QPSK.Ts = 1/RadPar.fs;                  % Symbol duration in seconds
QPSK.NumberofQPSK = 1;                  % Number of QPSK signals
QPSK.latShift = 0.01;                   % Shift of QPSK Transmitter from GRP latitude
QPSK.lonShift = -0.028;                   % Shift of QPSK Transmitter from GRP longitude
QPSK.Gain = 1;                          % Asumme omin-directional isotropic QPSK transmitter G = 1
QPSK.SIR = -20;                          % Assumed SIR in [dB]
%% Interfering Radar Signal Parameters
IR.fc = 0.1*RadPar.fo;                  % Carrier frequency same as SAR [Hz]
IR.ts = RadPar.ts;                      % [s] Sample time same as SAR
IR.bw = RadPar.bw;                      % [Hz] Bandwidth of the RF signal to calcualte the ramp rate
IR.T = RadPar.T;                        % [s] Pulse width
IR.K =  (IR.bw /IR.T);                  % [Hz/s] Ramp (chirp) rate
IR.t = -IR.T:IR.ts:IR.T;            % Time base vector for a single signal
% IR.t = FastTime;                        % Time base vector for signal
IR.NumberofIR = 2^5;                      % Number of Radar signals
IR.latShift = 0.02;                      % Shift of Radar Transmitter from GRP latitude
IR.lonShift = -0.01;                      % Shift of Radar Transmitter from GRP longitude
IR.Gain = 1;                            % Asumme omin-directional isotropic Radar transmitter G = 1
IR.SIR = 90;                            % Assumed SIR in [dB]