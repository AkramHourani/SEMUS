function [a] = FP04_GetGroundReflect(Targetlat,Targetlon,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2)
%
% To get the ground reflectivity based on the optical data 
%
% 1. Extract the edges of the swath region
% input : 
%   TargetLat = 11 x 11         | TargetLon = 11 x 11
%    [-33.7873 ..... -33.8356]      [ 151.3131 ..... 151.2428]
%   latSawthL1SoI = 1600 x 1    | lonSwathL1SoI =  1600 x 1   
%   latSawthL2SoI = 1600 x 1    | lonSwathL2SoI =  1600 x 1 
%
%   min([latSwathL1SoI;latSwathL2SoI]) = -33.8381
%   max([latSawthL1SoI;latSawthL2SoI]) = -33.7848 
%
%   min([lonSwathL1SoI;lonSwathL2SoI]) = 151.2428
%   max([lonSawthL1SoI;lonSawthL2SoI]) = 151.3131
% 
%   Hence, 
%   latEdge    = [-33.8381, -33.7848]
%   lonEdge    = [151.2428, 151.3131]
%
%      This is edge -->   o---o---o---o---o---o---o---o---o---o---o L1 
%                         |   |   |   |   |   |   |   |   |   |   |
%                         o   o   o   o   o   o   o   o   o   o   o
%                         .   .   .   .   .   .   .   .   .   .   .
%                         o   o   o   o   o   X   o   o   o   o   o
%                         .   .   .   .   .   .   .   .   .   .   .
%                         o   o   o   o   o   o   o   o   o   o   o
%                         |   |   |   |   |   |   |   |   |   |   |
%     This is edge -->    o---o---o---o---o---o---o---o---o---o---o L2
%
%       latEdge = [-33.8381; -33.7848]
%       lonEdge = [-33.8381; -33.7848]
%
%       
latEdge    = [min([latSawthL1;latSawthL2]), max([latSawthL1;latSawthL2])];
lonEdge    = [min([lonSwathL1;lonSwathL2]), max([lonSwathL1;lonSwathL2])];
FigOptical = figure('Position',[0 0 1000 1000]); % Create blank figure of 1000 pixels
clf
geoplot(Targetlat(:),Targetlon(:),'O');

%% Swath edge line

% Function : geolimits(latlim,lonlim) adjusts the limits of the current geographic axes, 
% geographic bubble chart, or map axes to include the latitude and longitude 
% limits specified by latlim and lonlim, respectively. 
% If the current axes is not a geographic axes, a geographic bubble chart, or a map axes, 
% or if there is no current axes, then the function creates a new geographic axes with the specified limits.
geolimits(latEdge,lonEdge)
geobasemap satellite
[latEdge, lonEdge] = geolimits; % Extract ground truth limits becuase the requested limits are not the same as the actual limits
%
% After geolimits : 
% latEdge = [-33.8394  -33.7835]
% lonEdge = [151.2334  151.3226] 
%
%
ax = gca;
ax.Position=[0 0 1 1];
%pause(5);
disp("If the image is loaded press any key to continue...")
pause
Optical = rgb2gray(frame2im(getframe(ax)));
close (FigOptical)
%
%   getframe(ax) = 
%   struct with fields:
%       cdata: [1436×2000×3 uint8]
%    colormap: []
%
%   frame2im(...) --> 1436×2000×3
%   frame2im : Return image data associated with movie frame
%
%   rgb2gray(...) --> 1436×2000×3
%   rgb2gray : Convert RGB image or colormap to grayscale
%
%% Interpolating the optical map
latVec = linspace(latEdge(1),latEdge(2),size(Optical,1));
%
% this will be : 
% latVec = linspace(-33.8394,-33.7835,1436) = [ -33.8394
% -33.8394 .... -33.7835] --> 1 x 1436 --> distribute into ... values
% lonVec = linspace(151.2334,151.3226,2000) = [151.2334	151.2334 ...
% 151.3226] --> 1 x 2000 --> distribute into ... values
lonVec = linspace(lonEdge(1),lonEdge(2),size(Optical,2));
%
% Create meshgrid lonvec --> row .... latvec --> column....
[lonVec,latVec] = meshgrid(lonVec,latVec);

% Function : B = flipud(A) returns A with its rows flipped in the up-down 
% direction (that is, about a horizontal axis).
latVec = flipud(latVec);

% Function : Vq = interp2(X,Y,V,Xq,Yq) returns interpolated values of 
% a function of two variables at specific query points (Xq, Yq) using linear interpolation. 
a = double(interp2(lonVec,latVec,Optical,Targetlon,Targetlat,'nearest'));

% this is a, 11 x 11
% 208 208 208 208	240	88	246	248	239	252	239
% 208 208 208 208	250	88	239	239	239	239	239
% 208 208 208 208	240	88	239	242	243	239	243
% 208 208 230 241	208	88	208	239	255	207	221
% 208 230 235 230	238	88	230	239	239	208	208
% 208 208 239 230	208	88	227	239	239	208	239
% 208 208 239 230	209	88	208	208	208	229	239
% 208 208 230 208	208	88	208	208	208	239	250
% 208 208 208 208	208	88	208	231	208	248	245
% 208 208 208 208	208	88	208	208	240	239	246
% 208 208 208 208	208	88	208	208	232	244	247


end
