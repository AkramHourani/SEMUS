%% Define radar signal
NI05a_RadarTxGain
% Generate basic chirp for the pulse duration
IRpulses = exp(1j*pi * IR.K * (IR.t-tauoIR).^2  )  .*(IR.t>(-IR.T/2) + tauoIR).*(IR.t<(IR.T/2) + tauoIR);
% IRpulses = exp(1j*pi * (-2 *IR.fc * IR.t + IR.K * IR.t.^2  ))  .*(IR.t>(-IR.T/2)).*(IR.t<(IR.T/2));
% imagesc(abs(IRpulses))
% Repeat this transmitted chirp for the duration of the SAR flight
RepeatFactor = round(time2num(Param.ScanDuration)/IR.T);
IRadarIQ = repmat(IRpulses.',1,RepeatFactor);
% spectrogram(Pulses_I)
% spectrogram(Pulses_I, 10, 0, 50, RadPar.fs, 'yaxis','centered');
% imagesc(real(IRpulses))
IRadarIQ = IRadarIQ.';
% imagesc(abs(IRadarIQ))
%% Define Timing
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all")); % Receiving window
TxTime = Param.tg-RxTime;                                           % Transmitting window= Pulse Repetition Interval (PRI)-Receiving window 

IRTime = length(IRadarIQ)/RadPar.fs;                                  % Radar signal transmission time
Discard = round(length(IRadarIQ)* TxTime/IRTime);                 % Number of samples on Radar data to be removed
Keep = round(length(IRadarIQ)* RxTime/IRTime);                    % Number of samples on Radar data to be considered
%% Creat radar Interference Matrix
% First define Numbers of section in radar data
sections = round(length(IRadarIQ) / (Discard + Keep));
% Creat Matrix
for i = 1: sections-1
  IRadarsignal(i,:) = IRadarIQ(i+(i*Discard)+((i-1)*Keep):i+(i*Discard)+(i*Keep));
end
% imagesc(real(IRadarsignal))
% IRSignal = IRadarsignal(1:size(sqd,1),1:size(sqd,2));
% speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.tg);       % Platform speed = sqrt(Param.mu/(h+Re))
Azimuthindex = linspace(-etaTotal*Param.tg*speed/2,etaTotal*Param.tg*speed/2,etaTotal);
[IR_Imgx, IR_Imgy, ~] = latlon2local(latIR,lonIR,0,GRP);

[~, IR_azimuth_index] = min(abs(Azimuthindex(:)-IR_Imgy));     % linear index of closest entry
IR_azimuth = ind2sub(size(Azimuthindex), IR_azimuth_index);   %// convert linear index to row and col

% Index of one end of the LoRa from the interferer tx location (LoRa_azimuth/2)
Index_before = IR_azimuth - floor(size(IRadarsignal,1)/2);
% Index of second end of the LoRa from the interferer tx location (LoRa_azimuth/2)
Index_after = IR_azimuth + floor(size(IRadarsignal,1)/2);

sIR1 = zeros(Index_before,size(sqd,2));
sIR2 = zeros(size(sqd,1)-Index_after,size(sqd,2));
sIRadar = cat(1,sIR1,IRadarsignal,sIR2);
IRSignal = sIRadar(1:size(sqd,1),:);
figure,imagesc(abs(IRSignal).^2)
sInfR = IRSignal .* sqrt(PIR);
%% Adjust Interference phase
% phase_IR = exp(-1i*2*pi*RadPar.Lambda/slantRangeIR);
% sInfR = sInfR .* repmat(phase_IR.',1,size(sInfR,2));
%% Plotting
% figure,imagesc(abs(sInfR).^2)
% rr= sInfR.';
% spectrogram(rr(:),500,0,500,RadPar.fs ,'yaxis','centered')
