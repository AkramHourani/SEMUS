%% Define radar signal gain wrt SARQPSK.fs
% Assume the radar transmitter location inside the swath
latIR = GRP(1) + IR.latShift;                             % Radar Tx latitude has same latitude as mid of swath
lonIR = GRP(2) + IR.lonShift;                             % Radar Tx longitude is out of the swath longitude by 0.2 degree
[azIR,elevIR,slantRangeIR] = geodetic2aer(latIR,lonIR,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRIR = elevIR + 90 - RadPar.AntOffNadir;
OffBoreSightAzIR = azIR - sataz;
% The zeta is added such that half the power is matching the beamwidth
zeta = 0.886;             
% Recieved power at SAR antenna from radar Tx
RxGainIR =single(RadPar.Gain * abs(sinc(OffBoreSightRIR*zeta/RadPar.BeamRange)).^2 ...
    .* abs(sinc(OffBoreSightAzIR*zeta/RadPar.BeamAz)).^2);
tauo_R = Ro/c;                                            % Delay of the Ground refernece point
tauoIR = single(slantRangeIR/c - tauo_R);                 % Single way delay from the radar
%% Received Power from the signal
PrdB = 10*log10(sum((abs(sqd)).^2,'all')/size(sqd,1));    % Received power [dB] using reference signal from all the dwell at GRP
% % For a required SIR
PIR_RxRef = PrdB - IR.SIR;                                % The desired reference value of IR Power in dBm - considered at the min distance
IR_DistRef = min(slantRangeIR);
% IRPower is adjustable  according to the desired SIR
IRPower = 10^(PIR_RxRef/20) .* fspl(IR_DistRef,RadPar.Lambda) ./ IR.Gain * max(sqrt(RxGainIR)); 
% % Convert to dB
IRPowerdB = 20*log10(IRPower);
% Received power at SAR radar from IR tx in watts
PIR = IRPower * IR.Gain * sqrt(RxGainIR)  ./ fspl(slantRangeIR,RadPar.Lambda);