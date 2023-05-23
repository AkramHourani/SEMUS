%% Antenna Dimensions
[~,~,Ro] = geodetic2aer(latSwathMid(MidEta),lonSwathMid(MidEta),0,Satlla(MidEta,1),Satlla(MidEta,2),Satlla(MidEta,3),E);   % Reference range at the centre of the swath at ground refernece point(GRP)
[~,~,R1] = geodetic2aer(latSwathL1(MidEta),lonSwathL1(MidEta),0,Satlla(MidEta,1),Satlla(MidEta,2),Satlla(MidEta,3),E);     % Reference range at near edge of the swath and the center of the dwell
[~,~,R2] = geodetic2aer(latSwathL2(MidEta),lonSwathL2(MidEta),0,Satlla(MidEta,1),Satlla(MidEta,2),Satlla(MidEta,3),E);     % Reference range at far edge of the swath and the center of the dwell
[~,~,R11] = geodetic2aer(latSwathL1(1),lonSwathL1(1),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);              % Reference range at near edge of the swath and the start of the dwell
[~,~,R12] = geodetic2aer(latSwathL1(end),lonSwathL1(end),0,Satlla(end,1),Satlla(end,2),Satlla(end,3),E);    % Reference range at near edge of the swath and end of the dwell
[~,~,R21] = geodetic2aer(latSwathL2(1),lonSwathL2(1),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);              % Reference range at far edge of the swath and the start of the dwell
[~,~,R22] = geodetic2aer(latSwathL2(end),lonSwathL2(end),0,Satlla(end,1),Satlla(end,2),Satlla(end,3),E);    % Reference range at far edge of the swath and the end of the dwell

[~,~,R1all] = geodetic2aer(latSwathLMid(:),lonSwathLMid(:),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);              % Reference range at near edge of the swath and the start of the dwell
grazang1a = grazingang(Param.h,R1all(1),'Curved');                     % Grazing angle (in degrees) from the scene center and the center of the dwell 
grazang1b = grazingang(Param.h,R1all(end),'Curved');                   % Grazing angle (in degrees) from the scene center and the center of the dwell 
incident1a = 90 - grazang1a;                                           % Incident angle of GRP
incident1b = 90 - grazang1b;                                           % Incident angle of GRP
PRF = 2 *Ro * abs((cos(incident1a)-cos(incident1b))) / RadPar.Lambda

% Doppler frequency calculation
speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.ts);
[sqang,~] = sarsquintang(SatECI,RadPar.AntOffNadir);                     % Squint angle for SAR data collection                        
Dopplerfre = 2 * speed * (RadPar.BeamAz * pi / 180) * cosd(max(sqang)) / RadPar.Lambda    % Doppler Frequency - from Jansing
Dopplerfreq = 4* max(DeltaR) / (RadPar.Lambda * time2num(Param.ScanDuration,"seconds"))   % PRF min - Doppler according to my calculation
PRF_min = 1.2 * Dopplerfreq
PRF_min = 1.2 * speed^2 * time2num(Param.ScanDuration,"seconds") / (RadPar.Lambda * Ro)
[~,~,R_nn_S] = geodetic2aer(latSwathL1(1),lonSwathL1(1),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);              % Reference range at near edge of the swath and the start of the dwell
[~,~,R_fn_S] = geodetic2aer(latSwathL2(1),lonSwathL2(1),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);              % Reference range at far edge of the swath and the start of the dwell
PRF_max = ((2 * (R_fn_S-R_nn_S) /c ) + RadPar.T )^ -1
[~,~,R_f_S] = geodetic2aer(latSwathL2(end),lonSwathL2(end),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);    % Reference range at far edge of the swath and the end of the dwell
[~,~,R_f_E] = geodetic2aer(latSwathL2(end),lonSwathL2(end),0,Satlla(end,1),Satlla(end,2),Satlla(end,3),E);    % Reference range at far edge of the swath and the end of the dwell
Dopplerfreq = 4* (R_f_S-R_f_E) / (RadPar.Lambda * time2num(Param.ScanDuration,"seconds"))   % PRF min - Doppler according to my calculation
PRF_min = 1.2 * Dopplerfreq

% System Parameters Calculation
% STEP.0 Angles calculation
grazango = grazingang(Param.h,Ro,'Curved');                               % Grazing angle (in degrees) from the scene center and the center of the dwell 
incidento = 90 - grazango                                           % Incident angle of GRP
grazang1 = grazingang(Param.h,R1,'Curved');                 % Grazing angle (in degrees) from the near end of the swath at the center of the dwell
incident1 = 90 - grazang1                                           % Near end incident angle
grazang2 = grazingang(Param.h,R2,'Curved');                 % Grazing angle (in degrees) from the far end of the swath at the center of the dwell
incident2 = 90 - grazang2                                           % Far end incident angle
theta_Emin = incident1 - (RadPar.AntOffNadir - RadPar.SwathWidthDeg/2)  % Earth angle corresponding to min incident angle
theta_Emax = incident2 - (RadPar.AntOffNadir + RadPar.SwathWidthDeg/2)  % Earth angle corresponding to max incident angle

% STEP.1 Range resolution calculation
slantrngres = bw2rangeres(RadPar.bw)                                % Convert bandwidth to Slant range resolution (m)
gndrngres =   slantrngres/sind(incident1)                           % Obtained Ground range resolution (m)

% STEP.2 Antenna dimensions
Ant_area = 4*speed*RadPar.Lambda*Ro*tand(RadPar.AntOffNadir)/c      % Antenna Area
Ant_height = 0.88 * RadPar.Lambda / (RadPar.BeamRange * pi /180 )    % Along height dimension / elevation - del -  which determine the swath width along Range
Ant_length = Ant_area / Ant_height
Ant_length = 0.88 * RadPar.Lambda / (RadPar.BeamAz * pi /180 )

% STEP.3 Azimuth resolution calculation
azres = Ant_length / 2                                              % Azimuth Resolution

% STEP.4 Swath width calculation
swathwidth = 2.5*Re * (theta_Emax - theta_Emin) * pi /180           % Ground swath width across Range direction
swathwidth = RadPar.Lambda * Ro *RadPar.SwathWidthDeg/ (Ant_height *RadPar.BeamRange/2*sind(grazango))  % Ground swath width across Range direction
swathwidth = (slantrange2(MidEta) - slantrange1(MidEta))/ sind((RadPar.AntOffNadir/2-RadPar.SwathWidthDeg))   % Ground swath width across Range direction
swathwidth = 4 * (slantrange2(MidEta) - slantrange1(MidEta))              % Ground swath width across Range direction==> prooved analytically
swathlength = speed * time2num(Param.ScanDuration)                                      % Ground swath length across Azimuth direction
% RadPar.Ns=1;
% swathlength = 2*RadPar.Ns*Ro*tand(RadPar.BeamAz/2)                  % Ground swath length across Azimuth direction if sub-swaths is selected

% STEP.5 PRF calculation
prf = sarprf(speed,Ant_length)
prfmax = c * Ant_height / (2 * Ro * RadPar.Lambda * tand(incidento) )
prfmin = 2 *speed  / Ant_length
[prfmin,prfmax] = sarprfbounds(speed,azres,swathlength,grazango,'PulseWidth',RadPar.T)



