%% Define LORA Message
N07_LoRaGain                                                        % Generate the LoRa gain and power according to the required SIR
message = "Hello World!" ;
%% Transmit Signal
signalIQ = N08_LoRaTx(message,LORA.BW,LORA.SF,LORAPowerdB,LORA.fs,LORA.Delta_f) ;
%% Define LORA Timing
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all")); % Receiving window
TxTime = Param.ts-RxTime;                                           % Transmitting window= Pulse Repetition Interval (PRI)-Receiving window 

LORATime = length(signalIQ)/fs;                                     % LoRa signal transmission time
Discard = round(length(signalIQ)* TxTime/LORATime);                 % Number of samples on LORA data to be removed
Keep = round(length(signalIQ)* RxTime/LORATime);                    % Number of samples on LORA data to be considered
%% Creat LORA Interference Matrix
% First define Numbers of section in LORA data
sections = round(length(signalIQ) / (Discard + Keep));
% Creat Matrix
for i = 1: sections-1
  slora(i,:) = signalIQ(i+(i*Discard)+((i-1)*Keep):i+(i*Discard)+(i*Keep));
end
% Adjust size of LORA signal to match the SAR raw data
sLORA = zeros(size(sqd,1),size(sqd,2));
for i = 1 : LORA.NumberofLoRa
scale(i) = 0.25 + (1 - 0.25) * rand(1,1);           % 0.25 is the min limit of the scale to ensure the data size is longer and 1 is the max limit of the scale -for SF =9, 0.8 for SF=12
LORA{i}(:,:) = padarray(slora,[(round(size(sqd,1)*scale(i))-size(slora,1)), (size(sqd,2)-size(slora,2))],0,'pre');
sLORAi = padarray(LORA{i}(:,:),[(size(sqd,1)-size(LORA{i}(:,:),1)), (size(sqd,2)-size(LORA{i}(:,:),2))],0,'post');
sLORA = sLORAi +sLORA;
end
sLORA = sLORA .* sqrt(PLORA);
% rr= sLORA.';
% spectrogram(rr(:),500,0,500,fs,'yaxis','centered')
% imagesc(real(sLORA))