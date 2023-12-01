% Assume the QPSK transmitter location inside the swath
latQPSK = GRP(1) + QPSK.latShift;                         % QPSK Tx latitude has same latitude as mid of swath
lonQPSK = GRP(2) + QPSK.lonShift;                         % QPSK Tx longitude is out of the swath longitude by 0.2 degree
[azQPSK,elevQPSK,slantRangeQPSK] = geodetic2aer(latQPSK,lonQPSK,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRQPSK = elevQPSK + 90 - RadPar.AntOffNadir;
OffBoreSightAzQPSK = azQPSK - sataz;

% Recieved power at SAR antenna from QPSK Tx
RxGainQPSK = single(RadPar.Gain * (sinc(OffBoreSightRQPSK*zeta/RadPar.BeamRange)).^2 ...
    .* (sinc(OffBoreSightAzQPSK*zeta/RadPar.BeamAz)).^2);

%% Received Power from the signal
PrdB = 20*log10(rms(sqd,2));                    % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PQPSK_RxRef = PrdB - QPSK.SIR;                            % The desired reference value of QPSK Power in dBm - considered at the min distance
QPSK_DistRef = min(slantRangeQPSK);
% QPSKPower is adjustable  according to the desired SIR
QPSKPower = 10^(PQPSK_RxRef/10) .* 10^(fspl(QPSK_DistRef,RadPar.Lambda)/10) ./ (QPSK.Gain * RxGainQPSK); 
% % Convert to dB
QPSKPowerdB = 10*log10(QPSKPower);
% Received power at SAR radar from QPSK tx in watts
PQPSK = QPSKPower * QPSK.Gain * RxGainQPSK  ./ 10^(fspl(slantRangeQPSK,RadPar.Lambda)/10);