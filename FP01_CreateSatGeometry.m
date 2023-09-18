function [SatECI,Satlla,DateVector] = FP01_CreateSatGeometry(startTime,stopTime,Param,Elem,satname)

% The satelliteScenario object represents a 3D arena consisting of satellites, 
% ground stations, and the interactions between them. 
% Use this object to model satellite constellations, 
% model ground station networks, 
% perform access analyses between the satellites and the ground stations, 
% and visualize the results.


%startTime = datetime('01-Jan-2022 08:00:00');   % [s] Set up the start time
%stopTime  = startTime + minutes(1) ; % [s] Set up the end time

sc  = satelliteScenario(startTime,stopTime,Param.ts);

% adds a Satellite object from Keplerian elements defined in the Geocentric Celestial Reference Frame (GCRF) to the satellite scenario.
sat = satellite(sc,Elem.a,Elem.e,Elem.Inc,Elem.RAAN,Elem.omega,Elem.TA,'OrbitPropagator','two-body-keplerian','name',satname);


% returns a 3-by-n-by-m array of the position history pos of each satellite in the vector sat, where n is the number of 
% time samples and m is the number of satellites. 
% The rows represent the x, y, and z coordinates of the satellite in the Geocentric Celestial Reference Frame (GCRF).
[SatECI,~,DateTime] = states(sat);

% converts the datetime or duration value t to a date vector—that is, 
% a numeric vector whose six elements represent the year, month, 
% day, hour, minute, and second components of t.
DateVector = datevec(DateTime);
DateVector(end,:)=[];       % Last sample is coming sometimes corrupted
SatECI(:,end)=[]; % Trim the last reading it has some errors. SatECI = Earth centered intertial

% converts Earth-centered inertial (ECI) coordinates, specified by position, to latitude, longitude, altitude (LLA) geodetic coordinates.
% The conversion is based on the Coordinated Universal Time (UTC) you specify.
Satlla = eci2lla(SatECI',DateVector);
% Satlla(:,3) = ones(size(Satlla,1),1)*Param.h; % Assume an ideal sphere, added 2023 to simplify the range migration process (it makes the approach distance to the center of the swath symetrical)


% % Visualize satellite 
% show(sat)
% groundTrack(sat,LeadTime=3600)
% %play(sc,PlaybackSpeedMultiplier=40)


%%uncomment to record video
% Create a Scenario and a Viewer
% viewer = satelliteScenarioViewer(sc,PlaybackSpeedMultiplier=40);
%
% % set camera
% camlat = Satlla(1,1);
% camlon = Satlla(1,2);
% camhgt = Satlla(1,3)*1.5;
% campos(viewer, camlat, camlon, camhgt)
% 
% % Specify the name of the video file
% videoFile = ['[0]'.satname '.mp4']; 
% % Create a VideoWriter object
% videoWriter = VideoWriter(videoFile, 'MPEG-4'); 
% % Set the frame rate of the video
% videoWriter.FrameRate = 30; 
% % Open the video writer
% open(videoWriter); 
% h = findall(0, 'Type', 'figure', 'Name', 'Satellite Scenario Viewer');
% for time = startTime:seconds(Param.ts):stopTime
%     FP091_Screencapture(h,'tempImg.png');
%     % Get the current frame as an image
%     frameImage = imread('tempImg.png'); 
%     % Write the frame to the video file
%     writeVideo(videoWriter, frameImage); 
%     % Advance to the next frame
%     viewer.CurrentTime=time; 
% end
% close(videoWriter);

close all force




end