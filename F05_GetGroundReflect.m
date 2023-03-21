function a = F05_GetGroundReflect(Targetlat,Targetlon,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2)
% Here we get the ground reflectivity based on the optical data 
% Extact the edges of the swath region
latEdge = [min([latSawthL1;latSawthL2]), max([latSawthL1;latSawthL2])];
lonEdge = [min([lonSwathL1;lonSwathL2]), max([lonSwathL1;lonSwathL2])];
FigOptical=figure('Position',[0 0 1000 1000]);                  % 1000 pixels
clf
% geoplot(Targetlat(:),Targetlon(:),'.');                        % Swath edge line - Only for demonstration, can not be used while processing
geolimits(latEdge,lonEdge)
geobasemap satellite
[latEdge, lonEdge] = geolimits;                                 % Extract ground truth limits becuase the requested limits are not the same as the actual limits
ax=gca;
ax.Position=[0 0 1 1];
disp("If the image is loaded press any key to continue...")
pause
Optical=rgb2gray(frame2im(getframe(ax)));
close (FigOptical)
%% Interpolating the optical map
latVec = linspace(latEdge(1),latEdge(2),size(Optical,1));
lonVec = linspace(lonEdge(1),lonEdge(2),size(Optical,2));
[lonVec,latVec]=meshgrid(lonVec,latVec);
latVec = flipud(latVec);
a = double(interp2(lonVec,latVec,Optical,Targetlon,Targetlat,'nearest'));
end