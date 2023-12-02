clc; clear; close all hidden
% This is the main file for simulating space SAR - Passive

%% 0. Load parameters
AP00_ParametersFerdi
Param.NtargetsAz = 11; % number of targets in each eta bin
Param.NtargetsRange = 11; % number of targets in each eta bin

%% 1. Passive SAR
%1a. Create Geometry setup (SoI / Satellite of Interest)
% This Script/function create the satellite orbit
[SatECISoI,SatllaSoI,DateVector]=   FP01_CreateSatGeometry(startTime,stopTime,Param,Elem,'Passive_SAR');
etaTotal                        =   length(DateVector); % Total number of slow time steps

%1b. Finding the swath
[latSwathMidSoI,lonSwathMidSoI,slantrangeMidSoI,SwathwidthSoI,latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI]=FP02_FindSwath(SatllaSoI,RadPar,E);

%% 2. Active SAR
%2a. Create Geometry setup (Interferer)
startTime = startTime + Delta_t;
stopTime  = startTime + Param.ScanDuration ;
%Elem.TA  = Elem.TA+0.015;  %  in degrees
%Elem.omega  = Elem.omega+0.015;
[SatECII,SatllaI,~] = FP01_CreateSatGeometry(startTime,stopTime,Param,Elem,'Active_SAR'); % This Script/function create the satellite orbit

%2b. Finding the swath of the SoI / Interferer
[latSwathMidI,lonSwathMidI,slantrangeMidI,Swathwidths_mI,latSwathL1I,lonSwathL1I,latSwathL2I,lonSwathL2I] = FP02_FindSwath(SatllaI,RadPar,E);
speed= mean(sqrt(sum((diff(SatECISoI,[],2)).^2)) /Param.ts);
swathlength = speed * time2num(Delta_t);                              % Ground swath length across Azimuth direction

%% 3. Find GRP
% This will find the GRP in the middle of the swath
% find the range migration of the middle of the swath
% This is the index of mid swath
MidEta = round(length(lonSwathL2SoI)/2);
[~,~,RSoI] = geodetic2aer(latSwathMidSoI(MidEta),lonSwathMidSoI(MidEta),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);
%RoSoI = min(RSoI);
GRP = [latSwathMidSoI(MidEta),lonSwathMidSoI(MidEta),0]; % ground reference point

% This is the distance between the Active SAR (Int) to the mid of the Passive SAR (SoI) swath
[~,~,RI] = geodetic2aer(latSwathMidSoI(MidEta),lonSwathMidSoI(MidEta),0,SatllaI(:,1),SatllaI(:,2),SatllaI(:,3),E);
%RoI = min(RI);
Ro = min(RSoI+RI);

%% Check the Doppler frequency by checking the maximum velocity of the swath corners
V_max = F10_VelocityCheck(latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI,SatllaSoI,E,Ro,Param);

%% Plotting
% Plot the Swath Center line of Passive SAR and Active SAR
% Plot the Swath Edge for Passive SAR and Active SAR
figure(1)
figure_1
% create_vid_1

%% Plot the Swath Edge for Passive SAR and Active SAR
figure(2)
figure_2
hold on

%% 4. Generate spatial sampling points (Targets)
% Note these are the targets with respect to the SoI (Passive SAR)
[Targetlat,Targetlon]= FP03_GenerateTargets(latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI,Param);

%% 5. Get ground reflectivity
a = FP04_GetGroundReflect(Targetlat,Targetlon,latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI);

figure(3)
figure_3
% create_vid_2

%% 6. Define The source satellite waveform
%
% [az,elev,slantRange] = geodetic2aer(lat,lon,h,lat0,lon0,h0,spheroid) 
% transforms the geodetic coordinates specified by lat, lon, and h 
% to the local azimuth-elevation-range (AER) spherical coordinates 
% specified by az, elev, and slantRange. 
% Specify the origin of the local AER system with the geodetic coordinates 
% lat0, lon0, and h0. Each coordinate input argument must match the others 
% in size or be scalar. Specify spheroid as the reference spheroid 
% for the geodetic coordinates.
%
%   midEta = 800 
%
%   [~,~,Edge1] = geodetic2aer( -33.8127, 151.3117,0, -33.8906,154.5327,5.1348e+05,E);
%   [~,~,Edge2] = geodetic2aer( -33.8102, 151.2443,0, -33.8906,154.5327,5.1348e+05,E);
%
%   Edge1 = 5.994e+05 ; Edge2 = 6.0279e+05 
%   Swathwidth_SARDistance = 3.3682e+03 m
%   SwathWidthTime = 2.2470e-05 s
%   FastTime = ( -1.3482e-05:1.6667e-08:1.3482e-05) --> 1 x 1618 
%   TimeLength = 1618
%   PulseWidthSamples = 300
%
[~,~,Edge1] = geodetic2aer(latSwathL1SoI(MidEta),lonSwathL1SoI(MidEta),0,SatllaSoI(MidEta,1),SatllaSoI(MidEta,2),SatllaSoI(MidEta,3),E);
[~,~,Edge2] = geodetic2aer(latSwathL2SoI(MidEta),lonSwathL2SoI(MidEta),0,SatllaSoI(MidEta,1),SatllaSoI(MidEta,2),SatllaSoI(MidEta,3),E);
% Then calculate Distance and Time
Swathwidth_SARDistance = abs(Edge1-Edge2);
SwathWidthTime = Swathwidth_SARDistance/c*2;
% Calculate FastTime
FastTime = (-SwathWidthTime/2*Param.Margin:RadPar.ts:SwathWidthTime/2*Param.Margin);
TimeLength = length(FastTime);
PulseWidthSamples = round(RadPar.T/(FastTime(end)-FastTime(1))*TimeLength);

%% 7. Generate base chirp (not necessary step, just for testing)
tau = 0;
sb  = exp(-1j*pi *   (2*RadPar.fo * tau - RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));%.*Window; % 1 x 1618

% Plotting the chirp
figure(4)
figure_4

%% 8. (Optional) you can select the Testing value for testing the script
% 0 = optical processing 
% 1 = GRP
% 2 = few targets testing
% 3 = unity reflection
%
% TargetLat --> 11 x 11
% TargetLon --> 11 x 11

Testing = 3 ; 
%FileName = 'SAR_Image1P.mat';

if Testing == 1 % this is for single targets testing
    Targetlat = GRP(1);
    Targetlon = GRP(2);
    a = 1;
    FileName = 'Test01P.mat';
end

NTesting = 5;
if Testing == 2 % this is for three targets testing
    ToPick = randsample(numel(Targetlat),NTesting) ;
    Targetlat = Targetlat(ToPick);
    Targetlon = Targetlon(ToPick);
    a = ones(NTesting,1);
    FileName = 'Test02P.mat';
end

if Testing == 3 % This will force the reflectivity to unity
    a = 1;
    FileName = 'Points.mat';
end

%% 9. Approx azimuth of the satellite to calculate the antenna pattern
%
%   satazSoI = azimuth(-33.8653,154.5331,-33.9160,154.5322,E) +90;
%   satazSoI = 270.7711
%
%   satazI = azimuth(-33.8653,154.5328,-33.9160,154.5320,E) +90;
%   satazI = 270.7711
%
if RadPar.Left == 1
    satazSoI = azimuth(SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(end,1),SatllaSoI(end,2),E) +90;
    satazI = azimuth(SatllaI(1,1),SatllaI(1,2),SatllaI(end,1),SatllaI(end,2),E) +90;
else
    satazSoI = azimuth(SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(end,1),SatllaSoI(end,2),E) -90;
    satazI = azimuth(SatllaI(1,1),SatllaI(1,2),SatllaI(end,1),SatllaI(end,2),E) -90;
end

%% 10. Reference sqd that will be used for template match filtering 
% tauoBi = (RoSoI+RoI)/c;% bi-static delay of the reference point
%
% Ro = 1.2022e+06 m
% tauoBi =  0.0040 s
% reference sqd is at GRP
%
tauoBi  = (Ro)/c;% bi-static delay of the reference point
sqd_ref = zeros(etaTotal,TimeLength);% initialize the reflection matrix, 1600 x 1618
sqd     = (zeros(etaTotal,TimeLength));
disp ('Generating the reference signal...')

for etaIdx=1:etaTotal % 1 : 1600
    sqd_ref(etaIdx,:) = FP05_CalcReflectionBi(1,GRP(1),GRP(2),SatllaI(etaIdx,:),SatllaSoI(etaIdx,:),RadPar,E,mean(satazI),mean(satazSoI),c,tauoBi,FastTime);
end

%% 11. Generate the reflected signal from the entire swath
% This is the longest part of the simulations
% the script will step through the azimuth (slow time) and generate the
% reflected signal from the entire swath
% sqd --> 1600 x 1618
tic
disp (['Starting simulation, total steps ',num2str(etaTotal)])
for etaIdx=1:etaTotal
    sqd(etaIdx,:) = FP05_CalcReflectionBi(a,Targetlat,Targetlon,SatllaI(etaIdx,:),SatllaSoI(etaIdx,:),RadPar,E,mean(satazI),mean(satazSoI),c,tauoBi,FastTime);
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

%% 12. Plot the raw unfocused SAR signal (Optional)
figure(5)
figure_5
% create_vid_3
save(FileName)
%save('Points.mat')
