%% Define the Message
NI02a_LoRaGain                                                      % Generate the transmitter gain and power according to the required SIR
message = char(randi([33 126],1,1e4));
%% Transmit Signal
LoRaIQ = NI02b_LoRaTx(message,LORA.BW,LORA.SF,LORAPowerdB,RadPar.fs,LORA.Delta_f) ;
% imagesc(abs(LoRaIQ).^2)
% spectrogram(LoRaIQ(:),500,0,500,RadPar.fs ,'yaxis','centered')
%% Define LORA Timing
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all")); % Receiving window of the platform
TxTime = Param.tg-RxTime;                                           % Transmitting window= Pulse Repetition Interval (PRI)-Receiving window 

LORATime = length(LoRaIQ)/RadPar.fs;                                % LoRa signal transmission time
Discard = round(length(LoRaIQ)* TxTime/LORATime);                 % Number of samples on LORA data to be removed
Keep = round(length(LoRaIQ)* RxTime/LORATime);                    % Number of samples on LORA data to be considered
%% Creat LORA Interference Matrix
% First define Numbers of section in LORA data
sections = round(length(LoRaIQ) / (Discard + Keep));
% Creat Matrix
for i = 1: sections-1
  slora(i,:) = LoRaIQ(i+(i*Discard)+((i-1)*Keep):i+(i*Discard)+(i*Keep));
end
sLORA1 = zeros(1500,size(sqd,2));
sLORA2 = zeros(345,size(sqd,2));
sLORA = cat(1,sLORA1,slora,sLORA2);

sLORA = sLORA .* sqrt(PLORA);
phase_LoRa = exp(-1i*2*pi*RadPar.Lambda/slantRangeLORA);
sLORA = sLORA .* repmat(phase_LoRa.',1,size(sLORA,2));
% rr= sLORA.';
% spectrogram(rr(:),500,0,500,RadPar.fs ,'yaxis','centered')
% imagesc(abs(sLORA).^2)