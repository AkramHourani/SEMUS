%% Find Recieved gain at SAR antenna from the RFI interferer location
% Assume the radar transmitter location inside the swath
latInterf = GRP(1) + IR.latShift;                             % Radar Tx latitude has same latitude as mid of swath
lonInterf = GRP(2) + IR.lonShift;                             % Radar Tx longitude is out of the swath longitude by 0.2 degree
[azInterf,elevInterf,slantRangeInterf] = geodetic2aer(latInterf,lonInterf,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRIR = elevInterf + 90 - RadPar.AntOffNadir;
OffBoreSightAzIR = azInterf - sataz;
% figure,plot(real(slantRangeIR))
% The zeta is added such that half the power is matching the beamwidth
zeta = 0.886;             
% Recieved gain at SAR antenna from radar Tx
RxGainIR =single(RadPar.Gain * abs(sinc(OffBoreSightRIR*zeta/RadPar.BeamRange)).^2 ...
                            .* abs(sinc(2*OffBoreSightAzIR*zeta/RadPar.BeamAz)).^2);
% figure;imagesc(abs(RxGainIR))   % figure;plot((10*log10(RxGainIR))) ; figure;plot(((RxGainIR)))
%% Received Power from the signal
PrdB = 20*log10(rms(sqd(:)));                             % Received power [dB] using reference signal from all the dwell at GRP
% figure,plot(PrdB);
% PrdB = 10*log10(max((abs(sqd(:)).^2)./length(sqd(:))));  % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR,
PrIRRxRef = PrdB - IR.SIR;                                % The desired reference received value of IR Power in dBm - considered at the min distance
% % Finding the Tx power from interferer based on the desired PIR_RXREF and SIR
% [IR_DistRef,indx] = min(slantRangeInterf);
% PrIRTxdB = PrIRRxRef + fspl(IR_DistRef,IR.Lambda) - 10*log10(IR.Gain) - 10*log10(RxGainIR(indx)); 
% % First find the average channel loss to maintain 
ChannelLossdB = -fspl(slantRangeInterf,IR.Lambda) + 10*log10(IR.Gain) + 10*log10(RxGainIR); 
ChannelLossdB_Av = 10*log10(mean(10.^(ChannelLossdB/10)));
PrIRTxdB = PrIRRxRef - ChannelLossdB_Av;
% Received power at SAR radar from this calculated IR Tx power
PrIRRxdB = PrIRTxdB + 10*log10(IR.Gain) + 10*log10(RxGainIR) - fspl(slantRangeInterf,RadPar.Lambda);
% Received power at SAR from the interferer in Watts
PrIRRx = 10.^(PrIRRxdB/10); % max(PrIRRxdB) 
% figure;imagesc(abs(PrIRRx))
PrIRRxdB_Av = 10*log10(mean(10.^(PrIRRxdB/10)));
PrSARdB = 20*log10(rms(sqd,2));
% figure,plot((PrIRRxdB));hold on;plot(PrSARdB);grid on;plot([1 etaTotal], PrdB *[1 1]);hold off
