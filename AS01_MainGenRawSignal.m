clc; clear; close all hidden 
% This is the main script for simulating space SAR
% The script will generate a raw SAR signal (baseband) based on the optial
% satellite image of the taregt swath
%% Load paratmers
% You can change paramters here
A00_Parameters
%% Create Geomtry setup
% This Scrip/function creat the satellite orbit
[SatECI] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem);
DateTime =  startTime:seconds(Param.dt):(startTime + seconds(Param.dt *(size(SatECI,2)-1)));
DateVector = datevec(DateTime);
GeoTime = 0:Param.dt:Param.dt *(size(SatECI,2)- 1);
%% Finding the swath
Satlla = eci2lla(SatECI',DateVector);
[latSawthMid,lonSwathMid,slantrangeMid,Swathwidths_m,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,sataz]=F02_FindSwath(Satlla,RadPar,E);
%% Find ground reference point
R=[];
for etaIdx=1:length(latSawthMid)
    [~,~,R(etaIdx,:)] = geodetic2aer(latSawthMid(etaIdx),lonSwathMid(etaIdx),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
    %plot(R(ctr,:))
    %hold on
end
% Find eta with most symetrical pass
[~,etaIdx]=min(abs(R(:,end)-R(:,1)));
R = R(etaIdx,:);
[Ro,Idx] = min(R);
GRP = [latSawthMid(Idx),lonSwathMid(Idx),0]; % ground reference point
%% Plot swath
geoplot(Satlla(:,1),Satlla(:,2)); %Satellite subline
hold on
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
SwathWidthTime = Swathwidths_m/c;
FastTime = (-SwathWidthTime/2:RadPar.ts:SwathWidthTime/2);
etaTotal=length(DateVector);
TimeLength = length(FastTime);
sqd=(zeros(etaTotal,TimeLength)); % initialize the reflection matrix
PulseWidthSamples = round(RadPar.T/(FastTime(end)-FastTime(1))*TimeLength);
% creating hanning window with zeros
%Window = [zeros(1,round((TimeLength-PulseWidthSamples)/2)),hann(round(PulseWidthSamples))'];
%Window = [Window,zeros(1,TimeLength-length(Window))];
%%   Generate base chrip (not nessasry step, just for testing)
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
Testing=0; % 0 for optical proccessing and 1 for 3 targets testing, and 2 for a single target
FileName = 'matlabOptical';
NTesting = 3;
if Testing==1 % this is for three targets testing
    ToPick =randsample(numel(Targetlat),NTesting) ;
    Targetlat = Targetlat(ToPick);
    Targetlon = Targetlon(ToPick);
    a = ones(NTesting,1);
end
if Testing==2 % this is for single targets testing
    Targetlat = latSawthMid(Idx);
    Targetlon = lonSwathMid(Idx);
    a = 1;
end
%% Reference sqd that will be used for template match filtering
disp ('Generating the reference signal...')
tauo = 2*Ro/c;% delay of the refernece point
parfor etaIdx=1:etaTotal
    sqd_ref(etaIdx,:) = F05_CalcReflection(a,GRP(1),GRP(2),Satlla(etaIdx,:),RadPar,E,sataz,c,tauo,FastTime);
end
%% This is the logest part of the simulations 
% the script will step through the azimuth (slow time) and generate the
% reflected signal from the entire swath
tic
disp (['Starting simulation, total steps ',num2str(etaTotal)])
parfor etaIdx=1:etaTotal
    sqd(etaIdx,:) =F05_CalcReflection(a,Targetlat,Targetlon,Satlla(etaIdx,:),RadPar,E,sataz,c,tauo,FastTime);
    disp(etaIdx)
end
toc
%% Plot the raw unfocused SAR signal (Optional)
figure(7)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(sqd));
pc.LineStyle='none';
colormap bone
xlabel('Fast time [ms]')
ylabel('Azimuth index')
title('Raw time domain (magnitude)')
%% Save the waveform
save('Test01.mat')
