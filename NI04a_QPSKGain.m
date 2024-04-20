% Assume the QPSK transmitter location inside the swath
latInterf = GRP(1) + QPSK.latShift;                         % QPSK Tx latitude has same latitude as mid of swath
lonInterf = GRP(2) + QPSK.lonShift;                         % QPSK Tx longitude is out of the swath longitude by 0.2 degree
[azInterf,elevInterf,slantRangeInterf] = geodetic2aer(latInterf,lonInterf,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRInterf = elevInterf + 90 - RadPar.AntOffNadir;
OffBoreSightAzInterf = azInterf - sataz;
if RadPar.Left == 0 % RadPar.Left == 0 for the case from South to North - RadPar.Left == 1 for the case from North to South 
    OffBoreSightAzIR = azInterf - sataz + 360;
    else
    OffBoreSightAzIR = azInterf - sataz;
end
% The zeta is added such that half the power is matching the beamwidth
zeta = 0.886;             
% Recieved power at SAR antenna from AM Tx
RxGainInterf = single(RadPar.Gain * (sinc(OffBoreSightRInterf*zeta/RadPar.BeamRange)).^2 ...
    .* (sinc(OffBoreSightAzInterf*zeta/RadPar.BeamAz)).^2);
% figure,imagesc(abs(RxGainInterf))
%% Received Power from the signal - Method 2
PrdB = 20*log10(rms(sqd(:)));           % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PSD_sqd = PrdB - 10*log10(RadPar.bw);   % Average PSD of SAR
PInterf_RxRef = PSD_sqd - (AM.SIR- 10*log10(AM.BW));  %The desired reference value of AM Power in dBm - considered at the min distance
% Interf_DistRef = min(slantRangeInterf);
%  First find the average channel loss to maintain average received power
ChannelLossdB = -fspl(slantRangeInterf,RadPar.Lambda) + 10*log10(AM.Gain) + 10*log10(RxGainInterf); 
ChannelLossdB_Av = 10*log10(mean(10.^(ChannelLossdB/10)));
PrInterfTxdB = PInterf_RxRef - ChannelLossdB_Av;
% Received power at SAR radar from this calculated IR Tx power
InterfPowerdB = PrInterfTxdB + 10*log10(AM.Gain) + 10*log10(RxGainInterf) - fspl(slantRangeInterf,RadPar.Lambda);
% Received power at SAR radar from AM tx in watts
PInterf = 10.^(InterfPowerdB/10); % RFI recieved power in watts
% figure,imagesc(abs(PInterf))
PrSARdB = 20*log10(rms(sqd,2));
% figure,plot((InterfPowerdB));hold on;plot(PrSARdB);grid on;plot([1 etaTotal], PrdB *[1 1]);hold off
%% Received Power from the signal
% PrdB = 20*log10(rms(sqd,2));                    % Received power [dB] using reference signal from all the dwell at GRP
% % % For a required SIR
% PQPSK_RxRef = PrdB - QPSK.SIR;                            % The desired reference value of QPSK Power in dBm - considered at the min distance
% QPSK_DistRef = min(slantRangeQPSK);
% % QPSKPower is adjustable  according to the desired SIR
% QPSKPower = 10^(PQPSK_RxRef/10) .* 10^(fspl(QPSK_DistRef,RadPar.Lambda)/10) ./ (QPSK.Gain * RxGainQPSK); 
% % % Convert to dB
% QPSKPowerdB = 10*log10(QPSKPower);
% % Received power at SAR radar from QPSK tx in watts
% PQPSK = QPSKPower * QPSK.Gain * RxGainQPSK  ./ 10^(fspl(slantRangeQPSK,RadPar.Lambda)/10);