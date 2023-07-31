%% Define LORA Message
NI02a_LoRaGain                                                      % Generate the LoRa gain and power according to the required SIR
message = "Hello World! Hello World! " ;
%% Transmit Signal
LoRaIQ = NI02b_LoRaTx(message,LORA.BW,LORA.SF,LORAPowerdB,LORA.fs,LORA.Delta_f) ;
%% Define LORA Timing
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all")); % Receiving window
TxTime = Param.tg-RxTime;                                           % Transmitting window= Pulse Repetition Interval (PRI)-Receiving window 

LORATime = length(LoRaIQ)/LORA.fs;                                % LoRa signal transmission time
Discard = round(length(LoRaIQ)* TxTime/LORATime);                 % Number of samples on LORA data to be removed
Keep = round(length(LoRaIQ)* RxTime/LORATime);                    % Number of samples on LORA data to be considered
%% Creat LORA Interference Matrix
% First define Numbers of section in LORA data
sections = round(length(LoRaIQ) / (Discard + Keep));
% Creat Matrix
for i = 1: sections-1
  slora(i,:) = LoRaIQ(i+(i*Discard)+((i-1)*Keep):i+(i*Discard)+(i*Keep));
end
% Adjust size of LORA signal to match the SAR raw data
sLORA = zeros(size(sqd,1),size(sqd,2));
for i = 1 : LORA.NumberofLoRa
    scale(i) = 0.25 + (1 - 0.25) * rand(1,1);           % 0.25 is the min limit of the scale to ensure the data size is longer and 1 is the max limit of the scale -for SF =9, 0.8 for SF=12
    LORAi{i}(:,:) = padarray(slora,[(round(size(sqd,1)*scale(i))-size(slora,1)); (size(sqd,2)-size(slora,2))],0,'pre');
    sLORAi = padarray(LORAi{i}(:,:),[(size(sqd,1)-size(LORAi{i}(:,:),1)); (size(sqd,2)-size(LORAi{i}(:,:),2))],0,'post');
    sLORA = sLORAi +sLORA;
end
sLORA = sLORA .* sqrt(PLORA);
phase_LoRa = exp(-1i*2*pi*RadPar.Lambda/slantRangeLORA);
sLORA = sLORA .* repmat(phase_LoRa.',1,size(sLORA,2));
% rr= sLORA.';
% spectrogram(rr(:),500,0,500,LORA.fs ,'yaxis','centered')
% imagesc(real(sLORA))