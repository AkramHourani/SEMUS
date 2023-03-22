function [Reflection] = F05_CalcReflection(a,Targetlat,Targetlon,Satlla,RadPar,E,sataz,c,tauo,FastTime)
% Obtain the geomtry for all points (targets) in the swath
[az,elev,slantRange] = geodetic2aer(Targetlat,Targetlon,0,Satlla(1),Satlla(2),Satlla(3),E);
OffBoreSightRange = elev+90 - RadPar.AntOffNadir;
OffBoreSightAz = az - sataz;
AntennaGain =single(RadPar.Gain * abs(sinc(OffBoreSightRange/RadPar.BeamRange*0.6)) .* abs(sinc(OffBoreSightAz/RadPar.BeamAz*0.6)));
tau = single(2*slantRange/c - tauo); % relative delay w. r. t. the middel of the swath (i.e to the reference point)

% Process refelctions from all points in the scene - Generate Quadrature
% Demodulated Signal Sqd - Time Domain - Equation 5.1 Book
Pulses = exp(1j*pi *   (-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2   )    ) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:)));
Reflection = sum(a(:).*AntennaGain(:)./single(slantRange(:)).^2.*Pulses,1);

end