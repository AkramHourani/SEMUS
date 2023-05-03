% function F14_RadarTx
%% Interfering Radar Signal Parameters
IR.fc = RadPar.fo;                      % Carrier frequency same as SAR [Hz]
IR.fs = RadPar.fs;                      % Sampling frequency same as SAR sampling frequency [Hz]
IR.ts = (1/IR.fs);                      % [s] Sample time
IR.bw = RadPar.bw;                      % [Hz] Bandwidth of the RF signal to calcualte the ramp rate
IR.T = RadPar.T;                        % [s] Pulse width
IR.K =  (IR.bw /IR.T);                  % [Hz/s] Ramp (chirp) rate
IR.t = FastTime;                        % Time base vector for the carrier modulated signal
% IR.t = -1e-2:1/IR.fs:1e-2;              % Time base vector for the carrier modulated signal
IR.NumberofIR = 1;                      % Number of Radar signals
IR.TxShift = 0.05;                      % Shift of Radar Transmitter from GRP longitude
IR.Gain = 1;                            % Asumme omin-directional isotropic Radar transmitter G = 1
IR.SIR = 20;                            % Assumed SIR in [dB]
%% Define radar signal
% Assume the radar transmitter location inside the swath
latIR = GRP(1);                                          % Radar Tx latitude has same latitude as mid of swath
lonIR = GRP(2) - IR.TxShift;                             % Radar Tx longitude is out of the swath longitude by 0.2 degree
[azIR,elevIR,slantRangeIR] = geodetic2aer(latIR,lonIR,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRIR = elevIR + 90 - RadPar.AntOffNadir;
OffBoreSightAzIR = azIR - sataz;
% The zeta is added such that half the power is matching the beamwidth
zeta = 50.76;             
AntennaGainIR =single(RadPar.Gain * abs(sinc(OffBoreSightRIR*pi/180*zeta/RadPar.BeamRange)).^2 .* abs(sinc(OffBoreSightAzIR*pi/180*zeta/RadPar.BeamAz)).^2);
tauo_R = Ro/c;                                          % Delay of the Ground refernece point
tauoIR = single(slantRangeIR/c - tauo_R);               % Single way delay from the radar
parfor idx=1:etaTotal
        Pulses_I = exp(1j*pi * (-2 *IR.fc * tauoIR(:) + IR.K*(IR.t-tauoIR(:)).^2  ) ) .*(IR.t>(-IR.T/2+tauoIR(:))).*(IR.t<(IR.T/2+tauoIR(:)));
        sqd_I(idx,:) = sum(sqrt(AntennaGainIR(:) .* sqrt(IR.Gain))./single(slantRangeIR(:)) .*Pulses_I,1);
end
%% Received Power from the signal
PrdB = 10*log10(sum((abs(sqd)).^2,'all')/size(sqd,1));   % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PIR_RxRef = PrdB - IR.SIR;                               % The desired reference value of IR Power in dBm - considered at the min distance
IR_DistRef = min(slantRangeIR);
% IRPower is adjustable  according to the desired SIR
IRPower = 10^(PIR_RxRef/20) .* fspl(IR_DistRef,RadPar.Lambda) ./ IR.Gain * max(sqrt(AntennaGainIR)); 
% % Convert to dB
IRPowerdB = 20*log10(IRPower);
% Received power at SAR radar from IR tx in watts
P_IR = IRPower * IR.Gain * sqrt(AntennaGainIR)  ./ fspl(slantRangeIR,RadPar.Lambda);%% Define Activation Timing
%% Total interference signal
sInfR = sqd_I .* sqrt(P_IR);
% imagesc(real(sIR))
% %% Define IR Timing
% TxTime = Param.ts;
% RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all"));
% 
% IRTime = length(sqd_I)/IR.fs;
% Discard = round(length(sqd_I)* TxTime/IRTime);         % Number of samples on IR data to be removed
% Keep= round(size(sqd_I,2)* RxTime/IRTime);                   % Number of samples on IR data to be considered
% %% Creat radar Interference Matrix
% % First define Numbers of section in radar data
% sections = round(length(sqd_I) / (Discard + Keep));
% % Creat Matrix
% for i = 1: sections-1
%   s_IR(i,:) = sqd_I(i+(i*Discard)+((i-1)*Keep):i+(i*Discard)+(i*Keep));
% end
% % Adjust size of radar signal to match the SAR raw data
% s_IR = zeros(size(sqd,1),size(sqd,2));
% for i = 1 : IR.NumberofIR
%     scale(i) = 0.25 + (1 - 0.25) * rand(1,1);           % 0.25 is the min limit of the scale to ensure the data size is longer and 1 is the max limit of the scale -for SF =9, 0.8 for SF=12
%     Radari{i}(:,:) = padarray(s_IR,[(round(size(sqd,1)*scale(i))-size(s_IR,1)); (size(sqd,2)-size(s_IR,2))],0,'pre');
%     sRadari = padarray(Radari{i}(:,:),[(size(sqd,1)-size(Radari{i}(:,:),1)); (size(sqd,2)-size(Radari{i}(:,:),2))],0,'post');
%     s_IR = sRadari +s_IR;
% end
% s_IR = s_IR .* sqrt(P_IR);
