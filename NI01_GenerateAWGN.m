%% Noise calculation
% % Method 1
% Pr = 10*log10(sum((abs(sqd_ref)).^2,'all')/size(sqd_ref,1));                    % Received power [dB] using reference signal from all the dwell at GRP

% % Method 2
% Pr = 10*log10(sum((abs(Power_ref)).^2,'all')/size(Power_ref,1));                % Received power [dB] using reference signal from mid the dwell at GRP

% Method 3
% STEP.1 Range resolution calculation
slantrngres = bw2rangeres(RadPar.bw);                                           % Convert bandwidth to Slant range resolution [m]
% STEP.2 Antenna dimensions
speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.tg);                       % Platform speed = sqrt(Param.mu/(h+Re))
Ant_length = RadPar.Lambda / (RadPar.BeamAz * pi /180 );                        % Antenna length / Azimuth dimension [m]
% STEP.3 Azimuth resolution calculation
azres = Ant_length / 2;                                                         % Azimuth Resolution [m]

Pr = 10*log10(RadPar.Pt * (10^(RadPar.Gain/10))^2 * RadPar.Lambda^2  ...        % Received power [dB] using radar equation
     * 10^(RadPar.NESO/10) * slantrngres * azres ./ ( (4*pi)^3 * Ro^4 ...
     * RadPar.kb * RadPar.AT *RadPar.Bn * 10^(RadPar.F/10) * 10^(RadPar.FeederL/10)));

% Generate noise N
P_N = Pr - Noise.SNR;                                                           % Noise power [dB]

% Generate the I/Q matrix of the noise 
P_N = 10.^(P_N/10);                                                             % Convert back to Watt
AWGN = sqrt(P_N)* (randn(size(sqd)) + 1i * (randn(size(sqd)) ) );               % Generate the AWGN the same size as recieved signal 