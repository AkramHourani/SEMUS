 clc; clear; close all hidden
% This is the main file for simulating space SAR
%% Load paratmers
A00_Parameters                                                                      % Script for parameters definitions
%% Create Geomtry setup - STEP1.Create Satellite Scenario
[SatECI,velocity,DateTime] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem);  % Script for creating the satellite scenario and orbit
DateVector = datevec(DateTime);                                                     % Convert datetime data into Date vector of 6 elements for the whole flight duration
GeoTime = 0:Param.ts:Param.ts *(size(SatECI,2)- 1);                                 % Geometrical sampling time - Azimuth sampling
%% Find the swath - STEP2.Geometric Simulator
[Satlla,latSawthMid,lonSwathMid,slantrangeMid,Swathwidths_m,latSawthL1,lonSwathL1,slantrange1,latSawthL2,lonSwathL2,slantrange2,sataz]=F02_FindSwath(SatECI,DateVector,RadPar,E);
%% Find the Ground Reference Point - GRP
[R,Ro,Idx,GRP] = F03_FindGRP(latSawthMid,lonSwathMid,Satlla,E);
%% Generate spatial sampling points (Tragets)
[Targetlat,Targetlon]= F04_GenerateTargets(latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,Param); % This is for optical-based targets
%% Get ground reflectrivity - STEP3.Reflectivity Simulator
F05_GetGroundReflect
% Converting to cartisian coordinates for plotting
[xEast,yNorth,~] = latlon2local(Targetlat,Targetlon,0,GRP);
%% Test antenna pattern - STEP4.Amplitude Simulator
[OffBoreSightRange, OffBoreSightAz] = meshgrid (-RadPar.BeamRange:0.001:RadPar.BeamRange,-RadPar.BeamAz:0.001:RadPar.BeamAz); 
% The zeta is added such that half the power is matching the beamwidth
zeta = 50.76;                                           % Empirically calculated
AntennaGain = RadPar.Gain * (sinc(OffBoreSightRange*pi/180*zeta/RadPar.BeamRange)).^2 .* (sinc(OffBoreSightAz*pi/180*zeta/RadPar.BeamAz)).^2;
%%  Generate the reference reflected waveform template s(eta,t)
% FastTime = -((slantrange2(Idx) - slantrange1(Idx))/c)*1.1:RadPar.ts:((slantrange2(Idx) - slantrange1(Idx))/c)*1.1;   % This will cover about 1.1*swath width
SwathWidthTime = Swathwidths_m/c;                                           % Time at mid of the swath using the distance
FastTime = (-SwathWidthTime/2)*1.1:RadPar.ts:(SwathWidthTime/2)*1.1;        % Range fasttime
SlowTime = 0:Param.ts:Param.ft;                                             % Azimuth slowtime
etaTotal=length(DateVector);                                                % Slowtime length
TimeLength = length(FastTime);                                              % Fasttime length
sqd=zeros(etaTotal,TimeLength);                                             % Initialize the reflection matrix
PulseWidthSamples = round(RadPar.T/(FastTime(end)-FastTime(1))*TimeLength);
%% creating hanning window with zeros
% Window = [zeros(1,round((TimeLength-PulseWidthSamples)/2)),hann(round(PulseWidthSamples))'];
% Window = [Window,zeros(1,TimeLength-length(Window))];
%%  Generate base chrip (not nessasry step, just for testing)
tau = 0;
% Hanning window added to the pulse
sb = exp(-1j*pi *   (2*RadPar.fo * tau - RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));%.*Window;
%% Select: Testing = 0 for optical / Testing = 1 for 3 targets / Testing = 2 for 1 target GRP
Testing= 0;                      % 0 for optical proccessing 

FileName = 'matlabOptical1';
NTesting = 3;

if Testing==1                   % 1 for three targets testing
    ToPick =randsample(numel(Targetlat),NTesting) ;
    Targetlat = Targetlat(ToPick);
    Targetlon = Targetlon(ToPick);
    a = ones(NTesting,1);
    FileName = 'matlabTestTargets';
end

if Testing==2                   % 2 for one target testing - GRP (Midswath point)
    Targetlat = latSawthMid(Idx);
    Targetlon = lonSwathMid(Idx);
    a = 1;
    FileName = 'matlabTestRef';
end
%% Genrate reflections s(eta,t) and Reference signal- STEP5.Waveform Generator
tauo = 2*Ro/c;                  % Delay of the GRP
% Reference sqd_ref that will be used for template match filtering
disp ('Generating the reference signal...')
parfor etaIdx=1:etaTotal
    sqd_ref(etaIdx,:) = F06_CalcReflection(a,latSawthMid(Idx),lonSwathMid(Idx),Satlla(etaIdx,:),RadPar,E,sataz,c,tauo,FastTime);
end

Power_ref = F06_CalcReflection(1,latSawthMid(Idx),lonSwathMid(Idx),Satlla(round(etaTotal/2),:),RadPar,E,sataz,c,tauo,FastTime);

% Scene reflections sqd - reflected signal from the entire swath
disp (['Starting simulation, total steps ',num2str(etaTotal)])
parfor etaIdx=1:etaTotal
    sqd(etaIdx,:) = F06_CalcReflection(a,Targetlat,Targetlon,Satlla(etaIdx,:),RadPar,E,sataz,c,tauo,FastTime);
    
%     figure(6)
%     plot(FastTime/1e-6,real(sqd(eta,:)))
%     xlabel('Time [\mus]')
%     ylabel('Reflection magnitude')
    disp(etaIdx)
%     waitbar(eta/etaTotal,Bar)
end
%%
close all hidden
save(FileName)

