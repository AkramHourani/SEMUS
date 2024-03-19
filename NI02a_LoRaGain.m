% Assume the LoRa transmitter location inside the swath
latLORA = GRP(1) + LORA.latShift;                         % LoRa Tx latitude shift w.r.t GRP
lonLORA = GRP(2) + LORA.lonShift;                         % LoRa Tx longitude shift w.r.t GRP
[azLORA,elevLORA,slantRangeLORA] = geodetic2aer(latLORA,lonLORA,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRLORA = elevLORA + 90 - RadPar.AntOffNadir;
OffBoreSightAzLORA = azLORA - sataz;
if RadPar.Left == 0 % RadPar.Left == 0 for the case from South to North - RadPar.Left == 1 for the case from North to South 
    OffBoreSightAzIR = azInterf - sataz + 360;
    else
    OffBoreSightAzIR = azInterf - sataz;
end
% The zeta is added such that half the power is matching the beamwidth
zeta = 0.886;             
% Recieved power at SAR antenna from LoRa Tx
RxGainLORA = single(RadPar.Gain * (sinc(OffBoreSightRLORA*zeta/RadPar.BeamRange)).^2 ...
    .* (sinc(OffBoreSightAzLORA*zeta/RadPar.BeamAz)).^2);
% figure,imagesc(abs(RxGainLORA))

%% Received Power from the signal
PrdB = 10*log10(sum((abs(sqd)).^2,'all')/size(sqd,1));                   % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PLORA_RxRef = PrdB - LORA.SIR;                                                % The desired reference value of LoRa Power in dBm - considered at the min distance
LoRa_DistRef = min(slantRangeLORA);
% LORAPower is adjustable  according to the desired SIR
LORAPower = 10^(PLORA_RxRef/20) .* fspl(LoRa_DistRef,RadPar.Lambda) ./ LORA.Gain * max(sqrt(RxGainLORA)); 
% % Convert to dB
LORAPowerdB = 20*log10(LORAPower);
% Received power at SAR radar from LORA tx in watts
PLORA = LORAPower * LORA.Gain * mean(sqrt(RxGainLORA))  ./ fspl(slantRangeLORA,RadPar.Lambda);
% figure,imagesc(abs(PLORA))