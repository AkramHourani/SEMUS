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
AM.fs = 2*RadPar.fo;                   % Sampling frequency same higher than SAR sampling frequency [Hz]
% AM.fs = 50*RadPar.fs;                   % Sampling frequency same higher than SAR sampling frequency [Hz]
AM.t = (0:1/AM.fs:0.2);                 % Time base vector for the carrier modulated signal
AM.fm = 1e6;                            % Baseband signal frequency [Hz]       
AM.NumberofAM = 1;                      % Number of AM signals
AM.TxShift = 0.05;                      % Shift of AM Transmitter from GRP longitude
AM.Gain = 1;                            % Asumme omin-directional isotropic AM transmitter G = 1
AM.SIR = 20;                            % Assumed SIR in [dB]
%% QPSK Signal Parameters
QPSK.NofBits=2^14;                      % Number of Transmitted bits
QPSK.fc = RadPar.fo;                    % Carrier frequency same as SAR sampling frequency [Hz]
QPSK.OF = 100;                          % Oversampling factor (multiples of fc) - at least 4 is better
QPSK.fs = RadPar.fs;                    % Sampling frequency same as SAR [Hz]
QPSK.NumberofQPSK = 1;                  % Number of QPSK signals
QPSK.TxShift = 0.05;                    % Shift of QPSK Transmitter from GRP longitude
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
LORA.TxShift = 0.05;                    % Shift of LoRa Transmitter from GRP longitude
LORA.Gain = 1;                          % Asumme omin-directional isotropic LORA transmitter G = 1
LORA.SIR = 20;                          % Assumed SIR in [dB]
