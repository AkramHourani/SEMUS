% Define AM Signal
NI03a_AMGain
% x = cos(2*pi*AM.t*AM.fm);
% signal_IQA = ((1+x) .*cos(2*pi*AM.fc*AM.t)).';
x = cos(2*pi*AM.t*AM.fm);
InterfIQ = (1+x).* (cos(2 * pi * AM.fc * AM.t)+ 1i *sin(2 * pi * AM.fc * AM.t));
% InterfIQ = ammod(x,AM.fc,AM.fs);                            % Double-sideband AM
% figure;imagesc(abs(InterfIQ))
% figure;plot(AM.t, InterfIQ);
%% Define acquisition window and Chopping the signal according to the acquisition window
HopStopWindow = round(seconds(Param.ScanDuration) / size(sqd,1) * RadPar.fs);  % The number of samples in the whole window
BasicFilter = [ones(1,size(sqd,2)) zeros(1, HopStopWindow-size(sqd,2))];
% figure;plot(BasicFilter)
Filterstream = logical(repmat(BasicFilter,1,size(sqd,1)));
% % Chopping the signal according to the acquisition window
% % Chopping the signal according to the acquisition window
if length(InterfIQ)>= length(Filterstream)
    InterfIQ = InterfIQ(1:length(Filterstream));
else
    InterfIQ = [InterfIQ,zeros(1, length(Filterstream) - length(InterfIQ))];  % Pad x with zeros at the end
end
% figure;plot(abs(IRstream(1:200000)));hold on; plot(Filterstream(1:200000))
% figure;imagesc(abs(LoRaIQ))
sInterf = InterfIQ(Filterstream);
% figure;plot(real(InterfSignal(1:200)));% figure;imagesc(abs(sLora))
sInterf = (reshape(sInterf,size(sqd.'))).';
% figure;imagesc(abs(sInterf))
%%  Apply the power
% figure,imagesc(abs(sLora).^2)
sAM = sInterf .* sqrt(PInterf);
% figure;imagesc(abs(sAM))















