% function F14_RadarTx
%% Interfering Radar Signal Parameters
IR.fc = RadPar.fo;                      % Carrier frequency same as SAR [Hz]
IR.fs = RadPar.fs;                      % Sampling frequency same higher than SAR sampling frequency [Hz]
IR.t = (0:1/IR.fs:0.2);                  % Time base vector for the carrier modulated signal
IR.fm = 1e6;                            % Baseband signal frequency [Hz]       
IR.NumberofAM = 1;                      % Number of Radar signals
IR.TxShift = 0.05;                      % Shift of Radar Transmitter from GRP longitude
IR.Gain = 1;                            % Asumme omin-directional isotropic Radar transmitter G = 1
IR.SIR = 20;                            % Assumed SIR in [dB]
%% Define radar signal
% Assume the radar transmitter location inside the swath
latIR = GRP(1);                                          % Radar Tx latitude has same latitude as mid of swath
lonIR = GRP(2) - IR.TxShift;                              % Radar Tx longitude is out of the swath longitude by 0.2 degree
[azR,elevR,slantRangeR] = geodetic2aer(latIR,lonIR,0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
OffBoreSightRIR = elevR + 90 - RadPar.AntOffNadir;
OffBoreSightAzIR = azR - sataz;
tauoR = (slantRangeR)/c;                                % Single way delay from the radar
parfor idx=1:etaTotal
    sqd(idx,:) =F05_CalcReflection(a,latIR,lonIR,Satlla(idx,:),RadPar,E,sataz,c,tauoR,FastTime);
end




%% Define Activation Timing
TxTime = Param.ts;
RxTime = abs(min(FastTime,[],"all")) + abs(max(FastTime,[],"all"));

AMTime = length(signal_IQA)/AM.fs;
Discard_Q = round(length(signal_IQA)* TxTime/AMTime);         % Number of samples on LORA data to be removed
Keep_Q = round(size(sqd,2)* RxTime/AMTime);                   % Number of samples on LORA data to be considered
