% Assume the LoRa transmitter location inside the swath
latLORA = GRP(1) + LORA.latShift;                         % LoRa Tx latitude shift w.r.t GRP
lonLORA = GRP(2) + LORA.lonShift;                         % LoRa Tx longitude shift w.r.t GRP
[azLORA,elevLORA,slantRangeLORA] = geodetic2aer(latLORA,lonLORA,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRLORA = elevLORA + 90 - RadPar.AntOffNadir;
OffBoreSightAzLORA = azLORA - sataz;
% The zeta is added such that half the power is matching the beamwidth
zeta = 0.886;             
% Recieved power at SAR antenna from LoRa Tx
RxGainLORA = single(RadPar.Gain * (sinc(OffBoreSightRLORA*zeta/RadPar.BeamRange)).^2 ...
    .* (sinc(OffBoreSightAzLORA*zeta/RadPar.BeamAz)).^2);
% figure,imagesc(abs(RxGainLORA))

%% Received Power from the signal
PrdB = 20*log10(rms(sqd,2));                    % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PLORA_RxRef = PrdB - LORA.SIR;                                                % The desired reference value of LoRa Power in dBm - considered at the min distance
LoRa_DistRef = min(slantRangeLORA);
% LORAPower is adjustable  according to the desired SIR
LORAPower = 10^(PLORA_RxRef/10) .* 10^(fspl(LoRa_DistRef,RadPar.Lambda)/10) ./ (LORA.Gain * RxGainLORA); 
% % Convert to dB
LORAPowerdB = 10*log10(LORAPower);
% Received power at SAR radar from LORA tx in watts
PLORA = LORAPower * LORA.Gain * RxGainLORA  ./ 10^(fspl(slantRangeLORA,RadPar.Lambda)/10);
% figure,imagesc(abs(PLORA))
