clc; clear; close all hidden 
% This is the main script for simulating space SAR
% The script will generate a raw SAR signal (baseband) based on the optial
% satellite image of the taregt swath
%% Load paratmers
% You can change paramters here
A00_Parameters
%% Create Geomtry setup - STEP1.Create Satellite Scenario
% This Scrip/function creat the satellite orbit
[SatECI,Satlla,DateVector] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem);
etaTotal=length(DateVector);                            % Total numeber of slow time steps
%% Finding the swath - STEP2.Geometric Simulator
[latSawthMid,lonSwathMid,slantrangeMid,Swathwidth,latSwathL1,lonSwathL1,latSwathL2,lonSwathL2]=F02_FindSwath(Satlla,RadPar,E);
%% This will find the GRP in the middle of the swath
%find the range migration of the middle of the swath
% This is the index of the mid of the swath across the dwell time
MidEta = round(length(lonSwathL2)/2);
% Find the reference range at the centre of the dwell at the ground refernece point (GRP)
[~,~,R] = geodetic2aer(latSawthMid(MidEta),lonSwathMid(MidEta),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
GRP = [latSawthMid(MidEta),lonSwathMid(MidEta),0];      % Ground Reference Point (GRP)
Ro = min(R);                                            % The reference range at the ground refernece point (GRP)
%% Plot swath
geoplot(Satlla(:,1),Satlla(:,2));                       % Satellite subline
hold on
 
geoplot(latSawthMid,lonSwathMid,'--');                  % Swath center line
geoplot(GRP(1),GRP(2),'x');                             % Swath center point
geoplot(latSwathL1,lonSwathL1,'color',ColorOrder(2,:)); % Swath edge line 1
geoplot(latSwathL2,lonSwathL2,'color',ColorOrder(2,:)); % Swath edge line 2
legend('satellite subtrack','swath mid track')
title('Swath location') 
drawnow 
%% Generate spatial sampling points (Tragets)
[Targetlat,Targetlon]= F03_GenerateTargets(latSwathL1,lonSwathL1,latSwathL2,lonSwathL2,Param); % This is for optical-based targets
%% Get ground reflectrivity - STEP3.Reflectivity Simulator
a = F04_GetGroundReflect(Targetlat,Targetlon,latSwathL1,lonSwathL1,latSwathL2,lonSwathL2);
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
%% Test antenna pattern (optional part of the script) - STEP4.Amplitude Simulator
figure(6)
[OffBoreSightRange, OffBoreSightAz] = meshgrid (-RadPar.BeamRange:0.1:RadPar.BeamRange,-RadPar.BeamAz:0.1:RadPar.BeamAz);
% The 1.2 is added such that half the power is matching the beamwidth
AntennaGain = RadPar.Gain * abs(sinc(OffBoreSightRange/RadPar.BeamRange*0.6)) .* abs(sinc(OffBoreSightAz/RadPar.BeamAz*0.6));
pc =pcolor(OffBoreSightAz,OffBoreSightRange,AntennaGain);
pc.LineStyle='none'; 
axis equal;
colorbar
xlabel('Azimuth direction [deg]')
ylabel('Range direction [deg]')
title('Antenna gain pattern example')
%%  Generate the reference reflected waveform template s(eta,t)
[~,~,Edge1] = geodetic2aer(latSwathL1(MidEta),lonSwathL1(MidEta),0,Satlla(MidEta,1),Satlla(MidEta,2),Satlla(MidEta,3),E);   % Range of the first edge of the swath
[~,~,Edge2]  = geodetic2aer(latSwathL2(MidEta),lonSwathL2(MidEta),0,Satlla(MidEta,1),Satlla(MidEta,2),Satlla(MidEta,3),E);  % Range of the second edge of the swath
Swathwidth_SARDistance = abs(Edge1-Edge2);                                                  % Swath width in meters
SwathWidthTime = Swathwidth_SARDistance/c*2;                                                % Swath time
FastTime = (-SwathWidthTime/2*Param.Margin:RadPar.ts:SwathWidthTime/2*Param.Margin);        % Range fasttime
TimeLength = length(FastTime);                                                              % Fasttime length
sqd=(zeros(etaTotal,TimeLength));                                                           % Initialize the reflection matrix
PulseWidthSamples = round(RadPar.T/(FastTime(end)-FastTime(1))*TimeLength);
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
Testing=0; % 0 for optical proccessing and 1 for GRP, 2 for few targets testing, and 3 for unity reflection
FileName = 'SAR_Image.mat';
if Testing==1           % This is for single targets testing
    Targetlat = GRP(1);
    Targetlon = GRP(2);
    a = 1;
    FileName = 'Test01.mat';

end

NTesting = 5;           % Defining number of testing targets
if Testing==2           % This is for Ntesting targets
    ToPick =randsample(numel(Targetlat),NTesting) ;
    Targetlat = Targetlat(ToPick);
    Targetlon = Targetlon(ToPick);
    a = ones(NTesting,1);
    FileName = 'Test02.mat';
end

if Testing==3            % This will force the reflectivity to unity
    a = 1;
    FileName = 'Test03.mat';
end
%% Approx azimuth of the satellite to clauclate the antenna pattern
if RadPar.Left == 1
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2),E) +90;
else
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2),E) -90;
end

%% Reference sqd_ref that will be used as template for matched filter
disp ('Generating the reference signal...')
tauo = 2*Ro/c;                              % Delay of the Ground refernece point
parfor eta=1:etaTotal
    [sqd_ref(eta,:)] = F05_CalcReflection(1,GRP(1),GRP(2),Satlla(eta,:),RadPar,E,sataz,c,tauo,FastTime);
end
%% This is the logest part of the simulations - STEP5.Waveform Generator
% Scene reflections sqd - reflected signal from the entire swath
% the script will step through the azimuth (slow time) and generate the reflected signal from the entire swath
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
save(FileName)
