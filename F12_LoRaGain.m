% Assume the LoRa transmitter location inside the swath
latLORA = GRP(1);                                        % LoRa Tx latitude has same latitude as mid of swath
lonLORA = GRP(2) - LORA.TxShift;                         % LoRa Tx longitude is out of the swath longitude by 0.2 degree
[azLORA,elevLORA,slantRangeLORA] = geodetic2aer(latLORA,lonLORA,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRLORA = elevLORA+90 - RadPar.AntOffNadir;
OffBoreSightAzLORA = azLORA - sataz;

zeta = 50.76;             

% Recieved power at SAR antenna from LoRa Tx
RxGainLORA = single(RadPar.Gain * (sinc(OffBoreSightRLORA*pi/180*zeta/RadPar.BeamRange)).^2 .* (sinc(OffBoreSightAzLORA*pi/180*zeta/RadPar.BeamAz)).^2);

% % Received Power from the signal
PrdB = 10*log10(sum((abs(sqd)).^2,'all')/size(sqd,1));                   % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PLORA_RxRef = PrdB - LORA.SIR;                                                % The desired reference value of LoRa Power in dBm - considered at the min distance
LoRa_DistRef = min(slantRangeLORA);
% LORAPower is adjustable  according to the desired SIR
LORAPower = 10^(PLORA_RxRef/20) .* fspl(LoRa_DistRef,RadPar.Lambda) ./ LORA.Gain * max(sqrt(RxGainLORA)); 
% % Convert to dB
LORAPowerdB = 20*log10(LORAPower);
% Received power at SAR radar from LORA tx in watts
PLORA = LORAPower * LORA.Gain * sqrt(RxGainLORA)  ./ fspl(slantRangeLORA,RadPar.Lambda);