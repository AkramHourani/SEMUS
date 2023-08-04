function [Reflection] = F06_CalcReflection(a,Targetlat,Targetlon,Satlla,RadPar,E,sataz,c,tauo,FastTime)
% Obtain the geomtry for all points (targets) in the swath
[az,elev,slantRange] = geodetic2aer(Targetlat,Targetlon,0,Satlla(1),Satlla(2),Satlla(3),E);
OffBoreSightRange = elev + 90 - RadPar.AntOffNadir;
OffBoreSightAz = az - sataz;
% The zeta is added such that half the power is matching the beamwidth
zeta = 0.886;             
% Use this code for GPU processing Mode
%============================================================================================================%
AntennaGain =gpuArray(single(RadPar.Gain * abs(sinc(OffBoreSightRange*zeta/RadPar.BeamRange)).^2 ...
    .* abs(sinc(OffBoreSightAz*zeta/RadPar.BeamAz)).^2));

% Relative delay w. r. t. the middle of the swath (i.e to the reference point)
tau = gpuArray(single(2*slantRange/c - tauo));      

% Process refelctions from all points in the scene - Generate Quadrature Demodulated Signal Sqd - Time Domain
Pulses = gpuArray(exp(1i*pi * (-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2  )  ) ...
    .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:))));
Reflection = gpuArray(sum(AntennaGain(:) .* a(:) .* sqrt(RadPar.Pt) ./single(slantRange(:)).^2 .* Pulses,1));
% Reflection = gpuArray(sum(AntennaGain(:).* RadPar.Lambda .* sqrt(RadPar.Pt) .* a(:) ...
%     ./ (single(slantRange(:)).^2 .* sqrt((4*pi)^3)) .* Pulses,1));

% Use this code for Parallel CPU processing Mode
%============================================================================================================%
% AntennaGain =single(RadPar.Gain * abs(sinc(OffBoreSightRange*zeta/RadPar.BeamRange)).^2 ...
    % .* abs(sinc(OffBoreSightAz*zeta/RadPar.BeamAz)).^2);

% Relative delay w. r. t. the middle of the swath (i.e to the reference point)
% tau = single(2*slantRange/c - tauo);      

% Process refelctions from all points in the scene - Generate Quadrature Demodulated Signal Sqd - Time Domain
% Pulses = exp(1i*pi *   (-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2   )    ) ...
%   .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:)));
% Reflection = sum(AntennaGain(:) .* a(:) ./single(slantRange(:)).^2 .* Pulses,1);
% Reflection = sum(AntennaGain(:).* RadPar.Lambda .* sqrt(RadPar.Pt) .* a(:) ...
%   ./ (single(slantRange(:)).^2 .* sqrt((4*pi)^3)) .* Pulses,1);
end