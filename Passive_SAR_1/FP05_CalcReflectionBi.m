function Reflection = FP05_CalcReflectionBi(a,Targetlat,Targetlon,SatllaSrc,SatllaDist,RadPar,E,satazSrc,satazDist,c,tauoBi,FastTime)

% Find the time advance (time to wait before capturing, this
% is because of the large difference between the mono-static range and the
% bi-static range of the GRP)
% TimeAdvance = tauoBi-tauoMono;

% this is a, 11 x 11
% 208 208 208 208	240	88	246	248	239	252	239
% 208 208 208 208	250	88	239	239	239	239	239
% 208 208 208 208	240	88	239	242	243	239	243
% 208 208 230 241	208	88	208	239	255	207	221
% 208 230 235 230	238	88	230	239	239	208	208
% 208 208 239 230	208	88	227	239	239	208	239
% 208 208 239 230	209	88	208	208	208	229	239
% 208 208 230 208	208	88	208	208	208	239	250
% 208 208 208 208	208	88	208	231	208	248	245
% 208 208 208 208	208	88	208	208	240	239	246
% 208 208 208 208	208	88	208	208	232	244	247
%
%   TargetLat = 11 x 11         | TargetLon = 11 x 11
%    [-33.7873 ..... -33.8356]      [151.3131 ..... 151.2428]
%
% FP05_CalcReflectionBi(a,Targetlat,Targetlon,SatllaI(etaIdx,:),SatllaSoI(etaIdx,:),RadPar,E,mean(satazI),mean(satazSoI),c,tauoBi,FastTime);
%
% SatllaSrc --> SatllaI(etaIdx,:) --> Active SAR [lat lon alt] (3 cols)
% SatllaDist --> SatllaSoI(etaIdx,:) --> Passive SAR [lat lon alt] (3 cols)
% satazSrc --> mean(satazI) --> 270.7711
% satazDist --> mean(satazSoI) --> 270.7711
% tauoBi -->  0.0040
% FastTime --> 1 x 1618 --> [-1.3482e-05 .... 1.3467e-05]
%
% RadPar.AntOffNadir = 30
%

%% 1. Source to target distance --> Active SAR to point target on the ground

% create_vid_FP05

% 1.a Find azimuth, elevation, and slantrange for each point target
[az,elev,slantRangeSrc] = geodetic2aer(Targetlat,Targetlon,0,SatllaSrc(1),SatllaSrc(2),SatllaSrc(3),E);

% slantRangeSrc = 601086.475505532
% elev = -60.0006916769104
% az = 270.235507102489 
% 
% Targetlat = -33.8114372387549
% Targetlon = 151.278094652670
%
% satllasc = [-33.8653350925734	154.532772197063	513467.275654927]
%

% 1.b Find Boresight in Range direction
OffBoreSightRange = elev + 90 - RadPar.AntOffNadir;
% -60.0006916769104 + 90 - 30 = -0.000691676910356875

% 1.c Find Boresight in Azimuth direction
OffBoreSightAz = az - satazSrc;
% 270.235507102489  - 270.771112629398 = -0.535605526908853

%
%                                |     OffBoreSightRange  pi        |2  |     OffBoreSightAz  pi        |2
% AntennaGainSrc = RadPar.Gain x |sinc(-----------------x----) x 0.6| x |sinc(--------------x----) x 0.6|
%   (Active SAR)                 |     RadPar.BeamRange   180       |   |     RadPar.BeamAz   180       |
%

% 1.d Finally find Antenna Gain
%AntennaGainSrc = gpuArray(single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180/RadPar.BeamRange*0.6)).^2 .* abs(sinc(OffBoreSightAz*pi/180/RadPar.BeamAz*0.6)).^2)); 
% 15.848866
AntennaGainSrc = single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180/RadPar.BeamRange*0.6)).^2 .* abs(sinc(OffBoreSightAz*pi/180/RadPar.BeamAz*0.6)).^2);

%% 2. Target to destination distance --> point target to Passive SAR

% 2.a Find azimuth, elevaiton, and slantrange for each point target
[az,elev,slantRangeDist] = geodetic2aer(Targetlat,Targetlon,0,SatllaDist(1),SatllaDist(2),SatllaDist(3),E);

% 2.b Find Boresight in Range direction
OffBoreSightRange = elev + 90 - RadPar.AntOffNadir;

% 2.c Find Boresight in Azimuth direction
OffBoreSightAz = az - satazDist;

%
%                                 |     OffBoreSightRange  pi        |2  |     OffBoreSightAz  pi        |2
% AntennaGainDist = RadPar.Gain x |sinc(-----------------x----) x 0.6| x |sinc(--------------x----) x 0.6|
%  (Passive SAR)                  |     RadPar.BeamRange   180       |   |     RadPar.BeamAz   180       |
%

% 2.d Finally find Antenna Gain
%AntennaGainDist =gpuArray(single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180/RadPar.BeamRange*0.6)).^2 .* abs(sinc(OffBoreSightAz*pi/180/RadPar.BeamAz*0.6)).^2));
AntennaGainDist =single(RadPar.Gain * abs(sinc(OffBoreSightRange*pi/180/RadPar.BeamRange*0.6)).^2 .* abs(sinc(OffBoreSightAz*pi/180/RadPar.BeamAz*0.6)).^2);

%% 3. Calculate delay
%
%              slantRangeSrc+slantRangeDist
% tau = single(---------------------------- - tauoBi)
%                           c
%
%tau = gpuArray(single((slantRangeSrc+slantRangeDist)/c - tauoBi)); % relative delay w. r. t. the middle of the swath (i.e to the reference point)
tau = single((slantRangeSrc+slantRangeDist)/c - tauoBi); % relative delay w. r. t. the middle of the swath (i.e to the reference point)
% tau = 1.0604032e-08 s

%% 4. Calculate reflections 
% Process reflections from all points in the scene 
% from equation 2.20 (fg)

%4.a Calculate Pulses
% Pulses = gpuArray(exp(1j*pi*(-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2)) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:))));
Pulses = exp(1j*pi*(-2*RadPar.fo * tau(:) + RadPar.K*(FastTime-tau(:)).^2)) .*(FastTime>(-RadPar.T/2+tau(:))).*(FastTime<(RadPar.T/2+tau(:)));
%                                                                
%                                              2               -T                   T
% Pulses = exp(-jpi(-2.fo.tau + K.(FastTime-tau))) x(FastTime > - +tau)x(FastTime < - +tau)
%                                                               2                   2
% 1 x 1618

%4.b Calculate Reflections
Reflection = sum(a(:).*AntennaGainSrc(:).*AntennaGainDist(:)./...
    (single(slantRangeSrc(:)).^2 + single(slantRangeDist(:)).^2).*Pulses,1);
%
%                       AntennaGainSrc x AntennaGainDist 
%  Reflection = sum(a x -------------------------------- x Pulses , 1)
%                                   2                  2
%                       slantRangeSrc  +  slantRangeDist 
%
 
% Reflection ---> sqd (1 x 1618) similar size to fasttime
end