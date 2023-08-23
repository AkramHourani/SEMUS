%% Define radar signal
NI05a_RadarTxGain
% Generate basic chirp for the pulse duration
IRpulses = exp(1j*pi * IR.K * IR.t.^2  )  .*(IR.t>(-IR.T/2)).*(IR.t<(IR.T/2));
% IRpulses = exp(1j*pi * (-2 *IR.fc * IR.t + IR.K * IR.t.^2  ))  .*(IR.t>(-IR.T/2)).*(IR.t<(IR.T/2));
% imagesc(real(IRpulses))
% Repeat this transmitted chirp for the duration of the SAR flight
RepeatFactor = round(time2num(Param.ScanDuration)/IR.T);
IRadarIQ = repmat(IRpulses,1,RepeatFactor);
% IRadarIQ = IRadarIQ.';
% imagesc(real(IRadarIQ))
%% Define LORA Timing
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
IRSignal = IRadarsignal(1:size(sqd,1),1:size(sqd,2));
% % Adjust size of radar signal to match the SAR raw data
% IRSignal = zeros(size(sqd,1),size(sqd,2));
% for i = 1 : IR.NumberofIR
%     scale(i) = 0.25 + (1 - 0.25) * rand(1,1);           % 0.25 is the min limit of the scale to ensure the data size is longer and 1 is the max limit of the scale -for SF =9, 0.8 for SF=12
%     IRi{i}(:,:) = padarray(IRadarsignal,[(round(size(sqd,1)*scale(i))-size(IRadarsignal,1)); (size(sqd,2)-size(IRadarsignal,2))],0,'pre');
%     sIR = padarray(IRi{i}(:,:),[(size(sqd,1)-size(IRi{i}(:,:),1)); (size(sqd,2)-size(IRi{i}(:,:),2))],0,'post');
%     IRSignal = sIR +IRSignal;
% end
sInfR = IRSignal .* sqrt(PIR);
% phase_IR = exp(-1i*2*pi*RadPar.Lambda/slantRangeLORA);
% sInfR = sInfR .* repmat(phase_IR.',1,size(sInfR,2));
%% Total interference signal
% sInfR = sqd_I .* sqrt(PIR);
% imagesc(abs(sInfR))