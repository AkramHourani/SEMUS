%% Antenna Dimensions
[~,~,Ro] = geodetic2aer(latSawthMid(Idx),lonSwathMid(Idx),0,Satlla(Idx,1),Satlla(Idx,2),Satlla(Idx,3),E);   % Reference range at the centre of the swath at ground refernece point(GRP)
[~,~,R1] = geodetic2aer(latSawthL1(Idx),lonSwathL1(Idx),0,Satlla(Idx,1),Satlla(Idx,2),Satlla(Idx,3),E);     % Reference range at near edge of the swath and the center of the dwell
[~,~,R2] = geodetic2aer(latSawthL2(Idx),lonSwathL2(Idx),0,Satlla(Idx,1),Satlla(Idx,2),Satlla(Idx,3),E);     % Reference range at far edge of the swath and the center of the dwell
[~,~,R11] = geodetic2aer(latSawthL1(1),lonSwathL1(1),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);              % Reference range at near edge of the swath and the start of the dwell
[~,~,R12] = geodetic2aer(latSawthL1(end),lonSwathL1(end),0,Satlla(end,1),Satlla(end,2),Satlla(end,3),E);    % Reference range at near edge of the swath and end of the dwell
[~,~,R21] = geodetic2aer(latSawthL2(1),lonSwathL2(1),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);              % Reference range at far edge of the swath and the start of the dwell
[~,~,R22] = geodetic2aer(latSawthL2(end),lonSwathL2(end),0,Satlla(end,1),Satlla(end,2),Satlla(end,3),E);    % Reference range at far edge of the swath and the end of the dwell

% Doppler frequency calculation
speed= mean(sqrt(sum((diff(pos,[],2)).^2)) /Param.ts);
[sqang,~] = sarsquintang(pos,RadPar.AntOffNadir);                     % Squint angle for SAR data collection                        
Dopplerfre = 2 * speed * (RadPar.BeamAz * pi / 180) * cosd(max(sqang)) / RadPar.Lambda % Doppler Frequency

% System Parameters Calculation
% STEP.0 Angles calculation
grazango = grazingang(h,Ro,'Curved');                               % Grazing angle (in degrees) from the scene center and the center of the dwell 
incidento = 90 - grazango                                           % Incident angle of GRP
grazang1 = grazingang(h,slantrange1(Idx),'Curved');                 % Grazing angle (in degrees) from the near end of the swath at the center of the dwell
incident1 = 90 - grazang1                                           % Near end incident angle
grazang2 = grazingang(h,slantrange2(Idx),'Curved');                 % Grazing angle (in degrees) from the far end of the swath at the center of the dwell
incident2 = 90 - grazang2                                           % Far end incident angle
theta_Emin = incident1 - (RadPar.AntOffNadir - RadPar.SwathBeam/2)  % Earth angle corresponding to min incident angle
theta_Emax = incident2 - (RadPar.AntOffNadir + RadPar.SwathBeam/2)  % Earth angle corresponding to max incident angle

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
swathwidth = RadPar.Lambda * Ro *RadPar.SwathBeam/ (Ant_height *RadPar.BeamRange/2*sind(grazango))  % Ground swath width across Range direction
swathwidth = (slantrange2(Idx) - slantrange1(Idx))/ sind((RadPar.AntOffNadir/2-RadPar.SwathBeam))   % Ground swath width across Range direction
swathwidth = 4 * (slantrange2(Idx) - slantrange1(Idx))              % Ground swath width across Range direction==> prooved analytically
swathlength = speed * Param.ft                                      % Ground swath length across Azimuth direction
swathlength = 2*RadPar.Ns*Ro*tand(RadPar.BeamAz/2)                  % Ground swath length across Azimuth direction if sub-swaths is selected

% STEP.5 PRF calculation
prf = sarprf(speed,Ant_length)
prfmax = c * Ant_height / (2 * Ro * RadPar.Lambda * tand(incidento) )
prfmin = 2 *speed  / Ant_length
[prfmin,prfmax] = sarprfbounds(speed,azres,swathlength,grazango,'PulseWidth',RadPar.T)



