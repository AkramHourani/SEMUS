% function F14_RadarTx
%% Interfering Radar Signal Parameters
IR.fc = RadPar.fo;                      % Carrier frequency same as SAR [Hz]
IR.fs = RadPar.fs;                      % Sampling frequency same as SAR sampling frequency [Hz]
IR.ts = (1/IR.fs);                      % [s] Sample time
IR.bw = 15e6;                               % [Hz] Bandwidth of the RF signal to calcualte the ramp rate
IR.T = 5e-6;                                % [s] Pulse width
IR.K =  (IR.bw /IR.T);              % [Hz/s] Ramp (chirp) rate
IR.t = (0:1/IR.fs:0.2);                 % Time base vector for the carrier modulated signal
IR.NumberofIR = 1;                      % Number of Radar signals
IR.TxShift = 0.05;                      % Shift of Radar Transmitter from GRP longitude
IR.Gain = 1;                            % Asumme omin-directional isotropic Radar transmitter G = 1
IR.SIR = 20;                            % Assumed SIR in [dB]
%% Define radar signal
% Assume the radar transmitter location inside the swath
latIR = GRP(1);                                          % Radar Tx latitude has same latitude as mid of swath
lonIR = GRP(2) - IR.TxShift;                              % Radar Tx longitude is out of the swath longitude by 0.2 degree
[azR,elevR,slantRangeR] = geodetic2aer(latIR,lonIR,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRIR = elevR + 90 - RadPar.AntOffNadir;
OffBoreSightAzIR = azR - sataz;
% The zeta is added such that half the power is matching the beamwidth
zeta = 50.76;             
AntennaGainR =single(IR.Gain * abs(sinc(OffBoreSightRIR*pi/180*zeta/RadPar.BeamRange)).^2 .* abs(sinc(OffBoreSightAzIR*pi/180*zeta/RadPar.BeamAz)).^2);

tauoR = (slantRangeR)/c;                                % Single way delay from the radar
parfor idx=1:etaTotal
        Pulses_I = exp(1j*pi * (-2 *IR.fc * tauoR(:) + IR.K*(IR.t-tauoR(:)).^2  ) ) .*(IR.t>(-IR.T/2+tauoR(:))).*(IR.t<(IR.T/2+tauoR(:)));
        sqd_I(idx,:) = sum(sqrt(AntennaGainR(:) * sqrt(IR.Gain))./single(slantRangeR(:)).^2.*Pulses_I,1);
end

%% Define Activation Timing
TxTime = Param.ts;
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all"));

IRTime = length(sqd_I)/IR.fs;
Discard_Q = round(length(sqd_I)* TxTime/IRTime);         % Number of samples on IR data to be removed
Keep_Q = round(size(sqd_I,2)* RxTime/IRTime);                   % Number of samples on IR data to be considered
