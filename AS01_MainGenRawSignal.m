clc; clear; close all hidden 
% This is the main script for simulating space SAR
% The script will generate a raw SAR signal (baseband) based on the optial
% satellite image of the taregt swath
%% Load paratmers
% You can change paramters here
A00_Parameters                                                                      % Script for parameters definitions
%% Create Geomtry setup - STEP1.Create Satellite Scenario
% This Scrip/function for creating the satellite scenario and orbit
[SatECI] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem);
DateTime =  startTime:seconds(Param.dt):(startTime + seconds(Param.dt *(size(SatECI,2)-1)));
DateVector = datevec(DateTime);                                                     % Convert datetime data into Date vector of 6 elements for the whole flight duration 
GeoTime = 0:Param.dt:Param.dt *(size(SatECI,2)- 2);                                 % Geometrical sampling time - Azimuth sampling
DateVector(end,:)=[];                                                               % Trim the last reading it has some errors
SatECI(:,end)=[];                                                                   % Trim the last reading it has some errors
%% Finding the swath
% Convert Earth-centered inertial (ECI) coordinates of satellite into latitude, longitude, altitude (LLA) geodetic coordinates
Satlla = eci2lla(SatECI',DateVector);                                               % The conversion is based on the Universal Coordinated Time (UTC) specified by Date vector
[latSawthMid,lonSwathMid,slantrangeMid,Swathwidth,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2]=F02_FindSwath(Satlla,RadPar,E);
%% Find the Ground Reference Point - GRP in the middle of the swath
% Find Index of mid point of the dwell
Idx = round(length(lonSwathL2)/2);
% Find the Reference range at the centre of the dwell at the GRP 
[~,~,R] = geodetic2aer(latSawthMid(Idx),lonSwathMid(Idx),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
[Ro,~] = min(R);
GRP = [latSawthMid(Idx),lonSwathMid(Idx),0];                                % Swath center point (Ground Reference Point GRP)

[~,~,SAR_Dist_Edge1] = geodetic2aer(latSawthL1(Idx),lonSwathL1(Idx),0,Satlla(Idx,1),Satlla(Idx,2),Satlla(Idx,3),E);
[~,~,SAR_Dist_Edge2] = geodetic2aer(latSawthL2(Idx),lonSwathL2(Idx),0,Satlla(Idx,1),Satlla(Idx,2),Satlla(Idx,3),E);
Swathwidth_SARDistance = abs(SAR_Dist_Edge1-SAR_Dist_Edge2);
%% Plot swath

geoplot(Satlla(:,1),Satlla(:,2)); %Satellite subline
hold on

geoplot(latSawthMid,lonSwathMid,'--'); %Swath center line

geoplot(latSawthMid,lonSwathMid,'--'); %Swath center line
geoplot(GRP(1),GRP(2),'x'); %Swath center point
geoplot(latSawthL1,lonSwathL1,'color',ColorOrder(2,:)); %Swath edge line 1
geoplot(latSawthL2,lonSwathL2,'color',ColorOrder(2,:)); %Swath edge line
legend('satellite subtrack','swath mid track')
title('Swath location')
drawnow
%% Generate spatial sampling points (Tragets)
[Targetlat,Targetlon]= F03_GenerateTargets(latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,Param); % This is for optical-based targets
%% Get ground reflectrivity
a = F04_GetGroundReflect(Targetlat,Targetlon,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2);
figure(2)
% Converting to cartisian coordinates for plotting
[xEast,yNorth,~] = latlon2local(Targetlat,Targetlon,0,GRP);
scatter(xEast(:)/1000,yNorth(:)/1000,2,a(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')
colormap bone
axis equal
hold on
plot(0,0,'+'); % Mid point (reference)
xlabel('x-axis [km]')
ylabel('y-axis [km]')
title('Satellite swath (optical)')
%% Test antenna pattern (optional part of the script)
figure(6)
[OffBoreSightRange,     OffBoreSightAz] = meshgrid (-5:5,-20:20);
% The 1.2 is added such that hal the power is matching the beamwidth
AntennaGain = RadPar.Gain * abs(sinc(OffBoreSightRange/RadPar.BeamRange*0.6)) .* abs(sinc(OffBoreSightAz/RadPar.BeamAz*0.6));
pc =pcolor(OffBoreSightAz,OffBoreSightRange,AntennaGain);
pc.LineStyle='none'; 
axis equal;
colorbar
xlabel('Azimuth direction [deg]')
ylabel('Range direction [deg]')
title('Antenna gain pattern example')
%%  Generate the reference reflected waveform template s(eta,t)
SwathWidthTime = Swathwidth_SARDistance/c*2;                                % Time at mid of the swath using the distance
FastTime = (-SwathWidthTime/2:RadPar.ts:SwathWidthTime/2);                  % Range fasttime
etaTotal=length(DateVector);                                                % Slowtime length
TimeLength = length(FastTime);                                              % Fasttime length
sqd=(zeros(etaTotal,TimeLength));                                           % Initialize the reflection matrix
PulseWidthSamples = round(RadPar.T/(FastTime(end)-FastTime(1))*TimeLength);
%% creating hanning window with zeros
%Window = [zeros(1,round((TimeLength-PulseWidthSamples)/2)),hann(round(PulseWidthSamples))'];
%Window = [Window,zeros(1,TimeLength-length(Window))];
%% Generate the base chrip (not nessasry step, just for testing)
tau = 0;
sb = exp(-1j*pi *   (2*RadPar.fo * tau - RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));%.*Window;
figure(5)
plot(FastTime/1e-6,real(sb))
xlabel('Time [\mus]')
ylabel('Real part')
title('reference pulse [mid swath point]')
drawnow
%% (Optional) you can select the Testing value for testing the script
% Select: Testing = 0 for optical / Testing = 1 for 3 targets / Testing = 2 for 1 target GRP
Testing=1;                                                                  % 0 for optical proccessing
FileName = 'matlabOptical';
NTesting = 3;

if Testing==1                                                               % (1) This is for three targets testing
    ToPick =randsample(numel(Targetlat),NTesting) ;
    Targetlat = Targetlat(ToPick);
    Targetlon = Targetlon(ToPick);
    a = ones(NTesting,1);
end

if Testing==2                                                               % (2) This is for single targets testing - GRP (Midswath point)
    Targetlat = GRP(1);
    Targetlon = GRP(2);
    a = 1;
end
%% Approx azimuth of the satellite to clauclate the antenna pattern
if RadPar.Left == 1
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2),E) +90;
else
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2),E) -90;
end
%% Genrate Reference signal sqd_ref(eta,t) and reflections sqd(eta,t) - STEP5.Waveform Generator
%% Reference sqd that will be used for template match filtering
disp ('Generating the reference signal...')
tauo = 2*Ro/c;                                                              % Delay of the GRP
parfor eta=1:etaTotal
    [sqd_ref(eta,:)] = F05_CalcReflection(a,GRP(1),GRP(2),Satlla(eta,:),RadPar,E,sataz,c,tauo,FastTime);
end
%% This is the longest part of the simulations 
% The script will step through the azimuth (slow time) and generate the
% reflected signal from the entire swath
tic
disp (['Starting simulation, total steps ',num2str(etaTotal)])
parfor eta=1:etaTotal
    sqd(eta,:) =F05_CalcReflection(a,Targetlat,Targetlon,Satlla(eta,:),RadPar,E,sataz,c,tauo,FastTime);
    disp(eta)
end
toc
%% Plot the raw unfocused SAR signal (Optional)
figure(7)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(sqd));
pc.LineStyle='none';
ax=gca;
grid on
ax.Layer = 'top';
colormap bone
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Raw time domain (magnitude)')
%% Save the waveform
save('Test02.mat')
