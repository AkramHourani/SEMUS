%% Noise and Interferers Parameters
%% AWGN signal parameter
Noise.SNR = 60;                         % Assumed Signal to noise ratio  [dB]
%% LoRa Signal Parameters
% Generate 0.6595s of LORA data ==> Controlled by the SF=11
% SF ∈ {7, 8, 9, 10, 11, 12}.           % Generate 0.16s of LORA data ==> Controlled by the SF=9
LORA.SF = 11;                           % Spreading factor
LORA.BW = 0.1e6 ;                       % Signal bandwidth of LoRa transmission [Hz]
% LORA.fc = RadPar.fo + 1e6 ;           % Carrier frequency [Hz]
LORA.Delta_f = 1e6 ;                    % Frequency Shift [Hz]
LORA.NumberofLoRa = 1;                  % Number of LoRa signals
LORA.latShift = 0.015;                   % Shift of LoRa Transmitter from GRP latitude
LORA.lonShift = -0.025;                 % Shift of LoRa Transmitter from GRP longitude
LORA.Gain = 1;                          % Asumme omin-directional isotropic LORA transmitter G = 1
LORA.SIR = 0;                          % Assumed SIR in [dB]
%% AM Signal Parameters
AM.fc = RadPar.fo;                      % Carrier frequency same as SAR [Hz]
AM.fs = 2*RadPar.fo;                    % Sampling frequency same higher than SAR sampling frequency [Hz]
AM.t = (0:1/AM.fs:0.2);                 % Time base vector for the carrier modulated signal
AM.fm = 1e6;                            % Baseband signal frequency [Hz]       
AM.NumberofAM = 1;                      % Number of AM signals
AM.latShift = 0.02;                     % Shift of AM Transmitter from GRP latitude
AM.lonShift = 0.02;                     % Shift of AM Transmitter from GRP longitude
AM.Gain = 1;                            % Asumme omin-directional isotropic AM transmitter G = 1
AM.SIR = 30;                            % Assumed SIR in [dB]
%% QPSK Signal Parameters
QPSK.NofBits=2^14;                      % Number of Transmitted bits
QPSK.fc = RadPar.fo;                    % Carrier frequency same as SAR sampling frequency [Hz]
QPSK.OF = 100;                          % Oversampling factor (multiples of fc) - at least 4 is better
QPSK.NumberofQPSK = 1;                  % Number of QPSK signals
QPSK.latShift = 0.02;                   % Shift of QPSK Transmitter from GRP latitude
QPSK.lonShift = 0.02;                   % Shift of QPSK Transmitter from GRP longitude
QPSK.Gain = 1;                          % Asumme omin-directional isotropic QPSK transmitter G = 1
QPSK.SIR = 30;                          % Assumed SIR in [dB]
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