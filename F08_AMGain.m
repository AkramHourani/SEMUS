% Assume the AM transmitter location inside the swath
latAM = GRP(1);                                          % AM Tx latitude has same latitude as mid of swath
lonAM = GRP(2) - AM.TxShift;                             % AM Tx longitude is out of the swath longitude by 0.2 degree
[azAM,elevAM,slantRangeAM] = geodetic2aer(latAM,lonAM,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRAM = elevAM + 90 - RadPar.AntOffNadir;
OffBoreSightAzAM = azAM - sataz;

% Recieved power at SAR antenna from AM Tx
RxGainAM = single(RadPar.Gain * (sinc(OffBoreSightRAM*pi/180*zeta/RadPar.BeamRange)).^2 .* (sinc(OffBoreSightAzAM*pi/180*zeta/RadPar.BeamAz)).^2);

% % Received Power from the signal
PrdB = 10*log10(sum((abs(sqd)).^2,'all')/size(sqd,1));         % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PAM_RxRef = PrdB - AM.SIR;                                     % The desired reference value of AM Power in dBm - considered at the min distance
AM_DistRef = min(slantRangeAM);
% AMPower is adjustable  according to the desired SIR
AMPower = 10^(PAM_RxRef/10) .* fspl(AM_DistRef,RadPar.Lambda) ./ AM.Gain * max(sqrt(RxGainAM)); 
% % Convert to dB
AMPowerdB = 10*log10(AMPower);
% Received power at SAR radar from AM tx in watts
PAM = AMPower * AM.Gain * sqrt(RxGainAM)  ./ fspl(slantRangeAM,RadPar.Lambda);