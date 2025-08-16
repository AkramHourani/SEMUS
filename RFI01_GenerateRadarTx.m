%% Step 1: Define radar signal
RFI02_RadarTxGain
% % Generate basic chirp for the pulse duration with Rect window
IRpulse = exp(1i * pi * IR.K * IR.t.^2  ) .*(IR.t > -IR.T/2).*(IR.t < IR.T/2);  
% % Generate basic chirp for the pulse duration with kaiser window
% IRpulse = exp(1j* pi * IR.K * IR.t.^2  ) .* (kaiser(length(IR.t),10)).';  
% Generate CW for simplicity
% IRpulse = exp(1j* 2 * pi * IR.freqShift * IR.t);IRpulses = IRpulse;   % ###         
% IRpulse = exp(1j* 2 * pi * IR.K * IR.t);        

% figure;spectrogram(IRpulse(:),4,0,256,IR.fs,'yaxis','centered');
% figure;plot(IR.t,real(IRpulse) );
% figure;pwelch(IRpulse,[],[],[],RadPar.fs,'centered')
%% Step 2: Apply the frequency shift by element-wise multiplication
IRpulses = IRpulse .*exp(1i * pi * 2 *  IR.freqShift * IR.t); % ###
% IRpulses = exp(1i * 2 * pi * IR.freqShift * IR.t); 
% Add AWGN to the chirp
% IRpulses = awgn(IRpulses,30); 
% figure;spectrogram(IRpulses(:),4,0,256,IR.fs ,'yaxis')
% figure;imagesc(IR.t,1:length(IRpulses),real(IRpulses))
% figure;pwelch(IRpulses);figure;pspectrum(IRpulses)
% figure;plot(IR.t,real((IRpulses)))
% figure;pwelch(IRpulses(:),[],[],[],RadPar.fs,'centered')
%% Step 3: Repeat this transmitted chirp for the duration of the flight
RepeatFactorGeom = ceil(seconds(Param.ScanDuration) * RadPar.fs/ length(IR.t));
% Repeat the chirp to ensure having enough samples before discarding
IRstream = repmat(IRpulses,1,RepeatFactorGeom);
% figure;imagesc(real(IRstream))
% figure;spectrogram(IRstream(:),4,0,256,RadPar.fs ,'yaxis')
% figure;pwelch(IRstream(:));figure;pspectrum(IRstream)
% figure;plot((abs((IRstream))))
%% Step 4: Define acquisition window and Chopping the signal according to the acquisition window
HopStopWindow = round(seconds(Param.ScanDuration) / size(sqd,1) * RadPar.fs);  % The number of samples in the whole window
BasicFilter = [ones(1,size(sqd,2)) zeros(1, HopStopWindow-size(sqd,2))];
% figure;plot(BasicFilter)
Filterstream = logical(repmat(BasicFilter,1,size(sqd,1)));
% figure;plot(Filterstream(1:200000))
% % Chopping the signal according to the acquisition window
IRstream = IRstream(1:length(Filterstream));
% figure;plot(abs(IRstream(1:200000)));hold on; plot(Filterstream(1:200000))
InterfSignal = IRstream(Filterstream);
% figure;plot(real(InterfSignal(1:200)));
InterfSignal = (reshape(InterfSignal,size(sqd.'))).';
%% Step 5: Apply a kaiser window over the rows of the RFI
% Window = kaiser(size(InterfSignal,2),1.5);
% RepeatWindow = repmat(Window.',size(InterfSignal,1),1);
% InterfSignalWindowed = InterfSignal .* RepeatWindow;
% figure;imagesc(abs(RepeatWindow));
% figure;imagesc(abs(InterfSignalWindowed))
%% Step 5: Define Doppler shift of the generated pulses
% % The SAR radial velocity wrt interferer
% slantRangeInterfI = (interp1(1:length(slantRangeInterf),slantRangeInterf,1:etaTotal+1,'pchip')).';
% V_IR = diff(slantRangeInterfI)/Param.tg;                  
% % Doppler Centroid frequency
% DopplerCentroid = - 2 * V_IR / IR.Lambda; 
% DopplerTime = linspace(0,2*max(FastTime),length(FastTime));
% for i= 1:etaTotal
%     DopplerSignal(i,:) = exp(1i * 2 * pi * DopplerCentroid(i) .* DopplerTime);
% end
% figure;imagesc(real(DopplerSignal))% figure;plot(real(DopplerSignal))
% figure;spectrogram(DopplerSignal(:),32,0,256,RadPar.fs,'yaxis','centered')% figure,pwelch(DopplerSignal(:))
%% Step 6: Apply the Radar Transmitter gain and power
sInterf =  InterfSignal .* sqrt(PrIRRx);
% figure;imagesc(abs(sqrt(PrIRRx)))
% figure;imagesc(abs(sInterf))
% figure;pwelch(sInterf(:))
% figure;imagesc(abs(sInterf))
% figure;plot((10*log10(sInterf)))
% sInterfT = sInterf.';
% figure;pwelch(InterfSignal(:),size(sqd,1),[],size(sqd,1),RadPar.fs,'centered')
