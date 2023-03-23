clc; clear; close all hidden 
% This is the main script for simulating space SAR
% The script will generate a raw SAR signal (baseband) based on the optial
% satellite image of the taregt swath
%% Load paratmers
% You can change paramters here
A00_Parameters
%% Create Geomtry setup
% This Scrip/function creat the satellite orbit
[SatECI,Satlla,DateVector] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem);
etaTotal=length(DateVector); % Total numeber of slow time steps
%% Finding the swath
[latSawthMid,lonSwathMid,slantrangeMid,Swathwidth,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2]=F02_FindSwath(Satlla,RadPar,E);
%% This will find the GRP in the middle of the swath
%find the range migration of the middle of the swath
% This is the index of mid swath
MidEta = round(length(lonSwathL2)/2);
[~,~,R] = geodetic2aer(latSawthMid(MidEta),lonSwathMid(MidEta),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
Ro = min(R);
GRP = [latSawthMid(MidEta),lonSwathMid(MidEta),0]; % ground reference point
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
[~,~,Edge1] = geodetic2aer(latSawthL1(MidEta),lonSwathL1(MidEta),0,Satlla(MidEta,1),Satlla(MidEta,2),Satlla(MidEta,3),E);
[~,~,Edge2]  = geodetic2aer(latSawthL2(MidEta),lonSwathL2(MidEta),0,Satlla(MidEta,1),Satlla(MidEta,2),Satlla(MidEta,3),E);
Swathwidth_SARDistance = abs(Edge1-Edge2);
SwathWidthTime = Swathwidth_SARDistance/c*2;
Margin = 1.2; % This is a margin to include targets farther than swatch width, i.e. for larger squinet angle targets
FastTime = (-SwathWidthTime/2*Margin:RadPar.ts:SwathWidthTime/2*Margin);
TimeLength = length(FastTime);
sqd=(zeros(etaTotal,TimeLength)); % initialize the reflection matrix
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
Testing=2; % 0 for optical proccessing and 1 for GRP, 2 for few targets testing, and 3 for unity reflection

if Testing==1 % this is for single targets testing
    Targetlat = GRP(1);
    Targetlon = GRP(2);
    a = 1;
    FileName = 'Test01.mat';

end

NTesting = 5;
if Testing==2 % this is for three targets testing
    ToPick =randsample(numel(Targetlat),NTesting) ;
    Targetlat = Targetlat(ToPick);
    Targetlon = Targetlon(ToPick);
    a = ones(NTesting,1);
    FileName = 'Test02.mat';
end

if Testing==3 % This will force the reflectivity to unity
    a = 1;
    FileName = 'Test03.mat';
end
%% Approx azimuth of the satellite to clauclate the antenna pattern
if RadPar.Left == 1
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2),E) +90;
else
    sataz = azimuth(Satlla(1,1),Satlla(1,2),Satlla(end,1),Satlla(end,2),E) -90;
end

%% Reference sqd that will be used for template match filtering
disp ('Generating the reference signal...')
tauo = 2*Ro/c;% delay of the refernece point
parfor eta=1:etaTotal
    [sqd_ref(eta,:)] = F05_CalcReflection(1,GRP(1),GRP(2),Satlla(eta,:),RadPar,E,sataz,c,tauo,FastTime);
end
%% This is the logest part of the simulations 
% the script will step through the azimuth (slow time) and generate the
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
save('Mesh_Mono.mat')
