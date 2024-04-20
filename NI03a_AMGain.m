% Assume the AM transmitter location inside the swath
latInterf = GRP(1) + AM.latShift;                             % AM Tx latitude has same latitude as mid of swath
lonInterf = GRP(2) + AM.lonShift;                             % AM Tx longitude is out of the swath longitude by 0.2 degree
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
%% Received Power from the signal -  Method 1
% PrdB = 10*log10(sum((abs(sqd)).^2,'all')/size(sqd,1));                   % Received power [dB] using reference signal from all the dwell at GRP
% % % For a required SIR
% PInterf_RxRef = PrdB - AM.SIR;                                                % The desired reference value of AM Power in dBm - considered at the min distance
% Interf_DistRef = min(slantRangeInterf);
% % AMPower is adjustable  according to the desired SIR
% InterfPower = 10^(PInterf_RxRef/20) .* fspl(Interf_DistRef,RadPar.Lambda) ./ AM.Gain * max(sqrt(RxGainInterf)); 
% % % Convert to dB
% PrInterfTxdB = 20*log10(InterfPower);
% % Received power at SAR radar from AM tx in watts
% PInterf = InterfPower * AM.Gain * mean(sqrt(RxGainInterf))  ./ fspl(slantRangeInterf,RadPar.Lambda);
% % figure,imagesc(abs(PInterf))