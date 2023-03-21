function Reflection = F06_CalcReflection(a,Targetlat,Targetlon,Satlla,RadPar,E,sataz,c,tauo,FastTime)

% Obtain the geometry for all points (targets) in the swath
[az,elev,slantRange] = geodetic2aer(Targetlat,Targetlon,0,Satlla(1),Satlla(2),Satlla(3),E);
OffBoreSightRange = elev + 90 - RadPar.AntOffNadir;
OffBoreSightAz = az - sataz;
% The zeta is added such that half the power is matching the beamwidth
zeta = 50.76;             
AntennaGain =single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180*zeta/RadPar.BeamRange)).^2 .* abs(sinc(OffBoreSightAz*pi/180*zeta/RadPar.BeamAz)).^2);

tau = single(2*slantRange/c - tauo);                % Relative delay w. r. t. the middel of the swath (i.e to the ground reference point)

%TimeLength = length(FastTime);
%PulseWidthSamples = round(RadPar.T/(FastTime(end)-FastTime(1))*TimeLength);
%Window = [zeros(1,round((TimeLength-PulseWidthSamples)/2)),hann(round(PulseWidthSamples))'];
%Window = [Window,zeros(1,TimeLength-length(Window))];
%Window = F06_ShiftArray(Window,-tau(:)/RadPar.ts);

% Process refelctions from all points in the scene - Generate Quadrature
% Demodulated Signal Sqd - Time Domain - Equation 5.1 Book
Pulses = exp(-1j*pi *   (2*RadPar.fo * tau(:) - RadPar.K*(FastTime-tau(:)).^2   )    ) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:)));
Reflection = sum(a(:).*sqrt(AntennaGain(:))./single(slantRange(:)).^2.*Pulses,1);
end