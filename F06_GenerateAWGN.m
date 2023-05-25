%% Noise calculation
% % Method 1
Pr = 10*log10(sum((abs(sqd_ref)).^2,'all')/size(sqd_ref,1));                    % Received power [dB] using reference signal from all the dwell at GRP

% % Method 2
% Pr = 10*log10(sum((abs(Power_ref)).^2,'all')/size(Power_ref,1));                % Received power [dB] using reference signal from mid the dwell at GRP

% % Method 3
% % STEP.1 Range resolution calculation
% slantrngres = bw2rangeres(RadPar.bw);                                           % Convert bandwidth to Slant range resolution [m]
% % STEP.2 Antenna dimensions
% speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.tg);                       % Platform speed = sqrt(Param.mu/(h+Re))
% Ant_area = 4*speed*RadPar.Lambda*Ro*tand(RadPar.AntOffNadir)/c;                 % Antenna Area [m2]
% Ant_height = 0.88 * RadPar.Lambda / (RadPar.BeamRange*(pi/180));                % Along height dimension / elevation - del -  which determine the swath width along Range [m]
% Ant_length = Ant_area / Ant_height;                                             % Antenna length / Azimuth dimension [m]
% % STEP.3 Azimuth resolution calculation
% azres = Ant_length / 2;                                                         % Azimuth Resolution [m]
% 
% Pr = 10*log10(Noise.Pt * (10^(RadPar.Gain/10))^2 * RadPar.Lambda^2  ...         % Received power [dB] using radar equation
%      * 10^(Noise.NESO/10) * slantrngres * azres ./ ( (4*pi)^3 * Ro^4 ...
%      * Noise.kb * Noise.T *Noise.Bn * 10^(Noise.F/10) * 10^(Noise.FeederL/10)));

% Generate noise N
P_N = Pr - Noise.SNR;                                                           % Noise power [dB]

% Generate the I/Q matrix of the noise 
P_N = 10.^(P_N/10);                                                             % Convert back to Watt
% AWGN = sqrt(P_N)* (randn(size(sqd)) + 1i * (randn(size(sqd)) ) );               % Generate the AWGN the same size as recieved signal 
AWGN = 0.01*sqrt(P_N)* (randn(size(sqd)) + 1i * (randn(size(sqd)) ) );          % Generate the AWGN the same size as recieved signal 
