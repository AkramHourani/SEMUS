function [Reflection] = F06_CalcReflection(a,Targetlat,Targetlon,Satlla,RadPar,E,sataz,c,tauo,FastTime)
% Obtain the geomtry for all points (targets) in the swath
[az,elev,slantRange] = geodetic2aer(Targetlat,Targetlon,0,Satlla(1),Satlla(2),Satlla(3),E);
OffBoreSightRange = elev + 90 - RadPar.AntOffNadir;
OffBoreSightAz = az - sataz;
% The zeta is added such that half the power is matching the beamwidth
zeta = 0.886;             
AntennaGain =gpuArray(single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180*zeta/RadPar.BeamRange)).^2 .* abs(sinc(OffBoreSightAz*pi/180*zeta/RadPar.BeamAz)).^2));

tau = gpuArray(single(2*slantRange/c - tauo));      % Relative delay w. r. t. the middel of the swath (i.e to the reference point)

% Process refelctions from all points in the scene - Generate Quadrature Demodulated Signal Sqd - Time Domain

% Use this code for GPU processing Mode
Pulses = gpuArray(exp(1j*pi *   (-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2   )    ) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:))));
Reflection = gpuArray(sum((a(:)) .* AntennaGain(:) ./single(slantRange(:)).^2 .* Pulses,1));
% Reflection = gpuArray(sum(AntennaGain(:).* RadPar.Lambda .* sqrt(RadPar.Pt) .* a(:) ./single(slantRange(:)).^2 .* sqrt((4*pi)^3) .* Pulses,1));

% % Use this code for Parallel CPU processing Mode
% Pulses = exp(1j*pi *   (-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2   )    ) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:)));
% Reflection = sum((a(:)) .* AntennaGain(:) ./single(slantRange(:)).^2 .* Pulses,1);
% Reflection = sum(AntennaGain(:).* RadPar.Lambda .* sqrt(RadPar.Pt) .* a(:) ./single(slantRange(:)).^2 .* sqrt((4*pi)^3) .* Pulses,1);
end