%% Define the Message  - Method 2
NI02a_LoRaGain                                                      % Generate the transmitter gain and power according to the required SIR
message = char(randi([22 126],1,1e4));
%% Transmit Signal
LoRaIQ = NI02b_LoRaTx(message,LORA.BW,LORA.SF,PrLORATxdB,RadPar.fs,LORA.Delta_f) ;
% figure;imagesc(abs(LoRaIQ))
% spectrogram(LoRaIQ(:),500,0,500,RadPar.fs ,'yaxis','centered')
LoRaIQ = LoRaIQ.';
%% Define acquisition window and Chopping the signal according to the acquisition window
HopStopWindow = round(seconds(Param.ScanDuration) / size(sqd,1) * RadPar.fs);  % The number of samples in the whole window
BasicFilter = [ones(1,size(sqd,2)) zeros(1, HopStopWindow-size(sqd,2))];
% figure;plot(BasicFilter)
Filterstream = logical(repmat(BasicFilter,1,size(sqd,1)));
% % Chopping the signal according to the acquisition window
% % Chopping the signal according to the acquisition window
if length(LoRaIQ)>= length(Filterstream)
    LoRaIQ = LoRaIQ(1:length(Filterstream));
else
    LoRaIQ = [LoRaIQ,zeros(1, length(Filterstream) - length(LoRaIQ))];  % Pad x with zeros at the end
end
% figure;plot(abs(IRstream(1:200000)));hold on; plot(Filterstream(1:200000))
% figure;imagesc(abs(LoRaIQ))
sLora = LoRaIQ(Filterstream);
% figure;plot(real(InterfSignal(1:200)));% figure;imagesc(abs(sLora))
sLora = (reshape(sLora,size(sqd.'))).';
% figure;imagesc(abs(sLora))
%%  Apply the power
% figure,imagesc(abs(sLora).^2)
sLORA = sLora .* sqrt(PLORA);
% figure;imagesc(abs(sLORA))
%% Adjust Interference phase
% phase_LoRa = exp(-1i*2*pi*RadPar.Lambda/slantRangeLORA);
% sLORA = sLORA .* repmat(phase_LoRa.',1,size(sLORA,2));
%% Plotting
% figure,imagesc(abs(sLORA).^2)
% rr= sLORA.';
% spectrogram(rr(:),500,0,500,RadPar.fs ,'yaxis','centered')
%% Define the Message  - Method 1
% NI02a_LoRaGain                                                      % Generate the transmitter gain and power according to the required SIR
% message = char(randi([33 126],1,1e4));
% %% Transmit Signal
% LoRaIQ = NI02b_LoRaTx(message,LORA.BW,LORA.SF,LORAPowerdB,RadPar.fs,LORA.Delta_f) ;
% % imagesc(abs(LoRaIQ).^2)
% % spectrogram(LoRaIQ(:),500,0,500,RadPar.fs ,'yaxis','centered')
% %% Define LORA Timing
% RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all")); % Receiving window of the platform
% TxTime = Param.tg-RxTime;                                           % Transmitting window= Pulse Repetition Interval (PRI)-Receiving window 
% 
% LORATime = length(LoRaIQ)/RadPar.fs;                                % LoRa signal transmission time
% Discard = round(length(LoRaIQ)* TxTime/LORATime);                 % Number of samples on LORA data to be removed
% Keep = round(length(LoRaIQ)* RxTime/LORATime);                    % Number of samples on LORA data to be considered
% %% Creat LORA Interference Matrix
% % First define Numbers of section in LORA data
% sections = round(length(LoRaIQ) / (Discard + Keep));
% % Creat Matrix
% for i = 1: sections-1
%   slora(i,:) = LoRaIQ(i+(i*Discard)+((i-1)*Keep):i+(i*Discard)+(i*Keep));
% end
% % speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.tg);       % Platform speed = sqrt(Param.mu/(h+Re))
% Azimuthindex = linspace(-etaTotal*Param.tg*speed/2,etaTotal*Param.tg*speed/2,etaTotal);
% [LORAImgx, LORAImgy, ~] = latlon2local(latLORA,lonLORA,0,GRP);
% 
% [~, LoRa_azimuth_index] = min(abs(Azimuthindex(:)-LORAImgy));     % linear index of closest entry
% LoRa_azimuth = ind2sub(size(Azimuthindex), LoRa_azimuth_index);   %// convert linear index to row and col
% 
% % Index of one end of the LoRa from the interferer tx location (LoRa_azimuth/2)
% Index_before = LoRa_azimuth - floor(size(slora,1)/2);
% % Index of second end of the LoRa from the interferer tx location (LoRa_azimuth/2)
% Index_after = LoRa_azimuth + floor(size(slora,1)/2);
% 
% sLORA1 = zeros(Index_before,size(sqd,2));
% sLORA2 = zeros(size(sqd,1)-Index_after,size(sqd,2));
% sLora = cat(1,sLORA1,slora,sLORA2);
% sLora = sLora(1:size(sqd,1),:);
% % figure,imagesc(abs(sLora).^2)
% sLORA = sLora .* sqrt(PLORA);
% %% Adjust Interference phase
% % phase_LoRa = exp(-1i*2*pi*RadPar.Lambda/slantRangeLORA);
% % sLORA = sLORA .* repmat(phase_LoRa.',1,size(sLORA,2));
% %% Plotting
% % figure,imagesc(abs(sLORA).^2)
% % rr= sLORA.';
% % spectrogram(rr(:),500,0,500,RadPar.fs ,'yaxis','centered')