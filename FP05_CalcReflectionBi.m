function Reflection = FP05_CalcReflectionBi(a,Targetlat,Targetlon,SatllaSrc,SatllaDist,RadPar,E,satazSrc,satazDist,c,tauoBi,FastTime)

% Find the time advance (time to wait before capturing, this
% is because of the large difference between the mono-static range and the
% bi-static range of the GRP)
% TimeAdvance=tauoBi-tauoMono;

% Source to target distance
[az,elev,slantRangeSrc] = geodetic2aer(Targetlat,Targetlon,0,SatllaSrc(1),SatllaSrc(2),SatllaSrc(3),E);

OffBoreSightRange = elev+90 - RadPar.AntOffNadir;

OffBoreSightAz = az - satazSrc;

AntennaGainSrc = gpuArray(single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180/RadPar.BeamRange*0.6)).^2 .* abs(sinc(OffBoreSightAz*pi/180/RadPar.BeamAz*0.6)).^2));
% AntennaGainSrc = single(RadPar.Gain * abs(sinc(OffBoreSightRange/RadPar.BeamRange*0.6)) .* abs(sinc(OffBoreSightAz/RadPar.BeamAz*0.6)));

% Target to destination distance
[az,elev,slantRangeDist] = geodetic2aer(Targetlat,Targetlon,0,SatllaDist(1),SatllaDist(2),SatllaDist(3),E);
OffBoreSightRange = elev+90 - RadPar.AntOffNadir;
OffBoreSightAz = az - satazDist;
AntennaGainDist =gpuArray(single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180/RadPar.BeamRange*0.6)).^2 .* abs(sinc(OffBoreSightAz*pi/180/RadPar.BeamAz*0.6)).^2));

tau = gpuArray(single((slantRangeSrc+slantRangeDist)/c - tauoBi)); % relative delay w. r. t. the middle of the swath (i.e to the reference point)
% tau = single((slantRangeSrc+slantRangeDist)/c - tauoBi); % relative delay w. r. t. the middle of the swath (i.e to the reference point)

% Process reflections from all points in the scene 
% from equation 2.20 (fg)
Pulses = gpuArray(exp(1j*pi*(-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2)) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:))));
% Pulses = exp(1j*pi*(-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2)) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:)));


Reflection = sum(a(:).*AntennaGainSrc(:).*AntennaGainDist(:)./...
    (single(slantRangeSrc(:)).^2 + single(slantRangeDist(:)).^2).*Pulses,1);
%%
end