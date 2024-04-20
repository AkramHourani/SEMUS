% Generate binary data stream
NI04a_QPSKGain
% Generate random bits
NofSymbols = QPSK.SymbolRate * time2num(Param.ScanDuration);
dataBits=randi([0 1],1,2 * QPSK.SymbolRate);   % Input binary data stream (0's and 1's) to modulate
% Reshape bits into pairs
dataSymbols = reshape(dataBits, 2, [])';
% Mapping of bits to QPSK symbols
% 00 -> 1+1i % 01 -> -1+1i % 11 -> -1-1i % 10 -> 1-1i
symbolMap = 100 * [1+1i, -1+1i, -1-1i, 1-1i];
mappedSymbols = symbolMap(bi2de(dataSymbols) + 1);
% % Generate baseband QPSK signal (Upsampling each symbol)
Ts = 1 / QPSK.SymbolRate; % Symbol duration in seconds
samplesPerSymbol = round(Ts * RadPar.fs);
qpskSignal = repmat(mappedSymbols, samplesPerSymbol, 1);
qpskSignal = qpskSignal(:).';
% Time vector for one symbol, assuming symbol rate is the same as sampling rate
t = (0:length(qpskSignal)-1) / RadPar.fs;
% Carrier signal
carrier = exp(1i*2*pi*QPSK.fc*t);
% Modulated signal
InterfIQ = real(qpskSignal .* carrier);
% figure;imagesc(abs(InterfIQ))
%% Define acquisition window and Chopping the signal according to the acquisition window
HopStopWindow = round(seconds(Param.ScanDuration) / size(sqd,1) * RadPar.fs);  % The number of samples in the whole window
BasicFilter = [ones(1,size(sqd,2)) zeros(1, HopStopWindow-size(sqd,2))];
% figure;plot(BasicFilter)
Filterstream = logical(repmat(BasicFilter,1,size(sqd,1)));
% % Chopping the signal according to the acquisition window
% % Chopping the signal according to the acquisition window
if length(InterfIQ)>= length(Filterstream)
    InterfIQnew = InterfIQ(1:length(Filterstream));
else
    InterfIQnew = [InterfIQ,zeros(1, length(Filterstream) - length(InterfIQ))];  % Pad x with zeros at the end
end
% figure;plot(abs(IRstream(1:200000)));hold on; plot(Filterstream(1:200000))
% figure;imagesc(abs(InterfIQnew))
sInterf = InterfIQnew(Filterstream);
% figure;plot(real(InterfSignal(1:200)));% figure;imagesc(abs(sLora))
sInterf = (reshape(sInterf,size(sqd.'))).';
% figure;imagesc(abs(sInterf))
%%  Apply the power
% figure,imagesc(abs(sLora).^2)
sQPSK = sInterf .* sqrt(PInterf);
% figure;imagesc(abs(sQPSK))
