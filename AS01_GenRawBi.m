clc; clear; close all hidden
% This is the main file for simulating space SAR - Passive
%% Load paratmers
addpath("SEMUS")
A00_Parameters
Param.NtargetsAz = 11; % number of targets in each eta bin
Param.NtargetsRange = 11; % number of targets in each eta bin

%% Create Geomtry setup (SoI)
% This Scrip/function creat the satellite orbit
[SatECISoI,SatllaSoI,DateVector] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem);
etaTotal=length(DateVector); % Total numeber of slow time steps
%% Finding the swath
[latSawthMidSoI,lonSwathMidSoI,slantrangeMidSoI,SwathwidthSoI,latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI]=F02_FindSwath(SatllaSoI,RadPar,E);
%% Create Geomtry setup (Interferer)
% Elem.Inc = 100;  % Inclination in degrees
% Elem.TA  = 124.735;  %  in degrees
% Elem.RAAN = 200.67;
% RadPar.T = 4e-6; % Pulse width
% RadPar.bw = 10e6; % Bandwidth of the RF signal to calcualte the ramp rate
% RadPar.K =  -(RadPar.bw /RadPar.T); % ramp rate
% RadPar.Left = 1; % the scanning on the left side of the satellite trajectory
Param.h = 300; 
[SatECII,SatllaI,~] = F01_CreateSatGeometry(startTime,stopTime,Param,Elem); % This Scrip/function creat the satellite orbit
%% Finding the swath of the SoI / Interferer
[latSawthMidI,lonSwathMidI,slantrangeMidI,Swathwidths_mI,latSwathL1I,lonSwathL1I,latSwathL2I,lonSwathL2I]=F02_FindSwath(SatllaI,RadPar,E);

%% This will find the GRP in the middle of the swath
%find the range migration of the middle of the swath
% This is the index of mid swath
MidEta = round(length(lonSwathL2SoI)/2);
[~,~,RSoI] = geodetic2aer(latSawthMidSoI(MidEta),lonSwathMidSoI(MidEta),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);
%RoSoI = min(RSoI);
GRP = [latSawthMidSoI(MidEta),lonSwathMidSoI(MidEta),0]; % ground reference point

% This is the distane between the interfering satellite to the mid of the
% SoI swath
[~,~,RI] = geodetic2aer(latSawthMidSoI(MidEta),lonSwathMidSoI(MidEta),0,SatllaI(:,1),SatllaI(:,2),SatllaI(:,3),E);
%RoI = min(RI);
Ro = min(RSoI+RI);
%% Plotting
Scale = 0.5;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);
geoplot(SatllaSoI(:,1),SatllaSoI(:,2),"color",ColorOrder(1,:),'LineWidth',2); %Satellite subline
hold on
geoplot(SatllaI(:,1),SatllaI(:,2),"color",ColorOrder(2,:),'LineWidth',2); %Satellite subline
geoplot(latSawthMidSoI,lonSwathMidSoI,':',"color",ColorOrder(1,:),'LineWidth',2); %Swath center line
geoplot(latSawthMidI,lonSwathMidI,':',"color",ColorOrder(2,:),'LineWidth',2); %Swath center line
geoplot(GRP(1),GRP(2),'x'); %Swath center point
geoplot(latSwathL1SoI,lonSwathL1SoI,'color',ColorOrder(1,:),'LineWidth',2); %Swath edge line 1
geoplot(latSwathL2SoI,lonSwathL2SoI,'color',ColorOrder(1,:),'LineWidth',2); %Swath edge line
geoplot(latSwathL1I,lonSwathL1I,'color',ColorOrder(2,:),'LineWidth',2); %Swath edge line 1
geoplot(latSwathL2I,lonSwathL2I,'color',ColorOrder(2,:),'LineWidth',2); %Swath edge line
geobasemap satellite
geolimits([-33.84  -33.74],[151.2 151.35])
set(gca,'fontsize',14,'Gridlinestyle','--');
set(gca,'LooseInset',get(gca,'TightInset'));
LG = legend('Satellite-of-interest','Interfering satellite','Fontsize',14);
LG.Position = [0.7    0.7375    0.2587    0.2126];
drawnow

%% Generate spatial sampling points (Tragets)
% Note these are the targets with respect to the SoI satellite
[Targetlat,Targetlon]= F03_GenerateTargets(latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI,Param);
%% Get ground reflectrivity
a = F04_GetGroundReflect(Targetlat,Targetlon,latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI);
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
title('Saetllite swath (optical)')
%% The source satellite waveform
[~,~,Edge1] = geodetic2aer(latSwathL1SoI(MidEta),lonSwathL1SoI(MidEta),0,SatllaSoI(MidEta,1),SatllaSoI(MidEta,2),SatllaSoI(MidEta,3),E);
[~,~,Edge2]  = geodetic2aer(latSwathL2SoI(MidEta),lonSwathL2SoI(MidEta),0,SatllaSoI(MidEta,1),SatllaSoI(MidEta,2),SatllaSoI(MidEta,3),E);
Swathwidth_SARDistance = abs(Edge1-Edge2);
SwathWidthTime = Swathwidth_SARDistance/c*2;
FastTime = (-SwathWidthTime/2*Param.Margin:RadPar.ts:SwathWidthTime/2*Param.Margin);
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
Testing=3; % 0 for optical proccessing and 1 for GRP, 2 for few targets testing, and 3 for unity reflection

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

    FileName = 'Mesh_Bi.mat';
end
%% Approx azimuth of the satellite to clauclate the antenna pattern
if RadPar.Left == 1
    satazSoI = azimuth(SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(end,1),SatllaSoI(end,2),E) +90;
    satazI = azimuth(SatllaI(1,1),SatllaI(1,2),SatllaI(end,1),SatllaI(end,2),E) +90;
else
    satazSoI = azimuth(SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(end,1),SatllaSoI(end,2),E) -90;
    satazI = azimuth(SatllaI(1,1),SatllaI(1,2),SatllaI(end,1),SatllaI(end,2),E) -90;
end
%% Reference sqd that will be used for template match filtering
%tauoBi = (RoSoI+RoI)/c;% bi-static delay of the refernece point
tauoBi = (Ro)/c;% bi-static delay of the refernece point
disp ('Generating the reference signal...')
parfor etaIdx=1:etaTotal
    sqd_ref(etaIdx,:) = F05_CalcReflectionBi(1,GRP(1),GRP(2),SatllaI(etaIdx,:),SatllaSoI(etaIdx,:),RadPar,E,mean(satazI),mean(satazSoI),c,tauoBi,FastTime);
end
%% This is the logest part of the simulations
% the script will step through the azimuth (slow time) and generate the
% reflected signal from the entire swath
tic
disp (['Starting simulation, total steps ',num2str(etaTotal)])
parfor etaIdx=1:etaTotal
    sqd(etaIdx,:) =F05_CalcReflectionBi(a,Targetlat,Targetlon,SatllaI(etaIdx,:),SatllaSoI(etaIdx,:),RadPar,E,mean(satazI),mean(satazSoI),c,tauoBi,FastTime);
    disp(etaIdx)
end
toc
%%
% close all
% for ctr=1:3
% [~,~,RTargetSoI] = geodetic2aer(Targetlat(ctr),Targetlon(ctr),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);
% [~,~,RTargetI] = geodetic2aer(Targetlat(ctr),Targetlon(ctr),0,SatllaI(:,1),SatllaI(:,2),SatllaI(:,3),E);
% ttt = (RTargetI+RTargetSoI)/c-tauoBi;
% plot(ttt)
% hold on
% end

%% Plot the raw unfocused SAR signal (Optional)
figure(7)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(sqd));
pc.LineStyle='none';
colormap bone
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Raw time domain (magnitude)')

%%
%save(FileName)
save('Mesh_Bi')
