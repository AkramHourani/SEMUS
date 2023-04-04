% Define AM Signal
% x = cos(2*pi*AM.t*AM.fm);
% signal_IQA = ((1+x) .*cos(2*pi*AM.fc*AM.t)).';
x = sin(2*pi*AM.t*AM.fm);
signal_IQA = ammod(x,AM.fc,AM.fs);                            % Double-sideband AM
%% Define AM Timing
TxTime = Param.ts;
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all"));

AMTime = length(signal_IQA)/AM.fs;
Discard_Q = round(length(signal_IQA)* TxTime/AMTime);         % Number of samples on LORA data to be removed
Keep_Q = round(size(sqd,2)* RxTime/AMTime);                   % Number of samples on LORA data to be considered
%% Creat AM Interference Matrix
% First define Number of section in AM data
sections = round(length(signal_IQA) / (Discard_Q + Keep_Q));
% Creat Matrix
for i = 1: sections-1
  sam(i,:) = signal_IQA(i+(i*Discard_Q)+((i-1)*Keep_Q):i+(i*Discard_Q)+(i*Keep_Q));
end
% Adjust size of AM signal to match the SAR raw data
sAM = zeros(size(sqd,1),size(sqd,2));
for i = 1 : AM.NumberofAM
scale(i) = 0.25 + (1 - 0.25) * rand(1,1);                     % 0.25 is the min limit of the scale and 1 is the max limit of the scale
AMi{i}(:,:) = padarray(sam,[(round(size(sqd,1)*scale(i))-size(sam,1)), 0],0,'pre');
sAMi = padarray(AMi{i}(:,:),[(size(sqd,1)-size(AMi{i}(:,:),1)), 0],0,'post');
sAM = sAMi +sAM;
end
% imagesc(real(sAM))
F08_AMGain
sAM = sAM .* sqrt(PAM);
