function [gh] = plot_Moon(npts,pos)
% Plots the Moon in the at the position in [km] as provided.
% 
% 
% Note: This function requires an image of the Moon (a texturemap) to plot on
% the ellipsoidal surface. You can get one you like from:
% https://svs.gsfc.nasa.gov/cgi-bin/details.cgi?aid=4720
% 
% The one included with this function is the map from July:
% https://svs.gsfc.nasa.gov/cgi-bin/details.cgi?aid=4720
% 
% Inputs:
%   npts: [](1x1)<integer> number of equatorial points to use to generate the
%         surface of the Moon
%   pos:  [](3,1)<double> cartesian position vector to center the moon plot
%         about
% 
% Outputs:
%   gh: [](1x1)<graphics handle> A graphics handle that contains data about the
%        plot if needed
% Equatorial and polar radii of the Moon:
RM_eq = 1738.1; %[km] Earth equatorial radius (NASA Fact Sheet)
RM_pol = 1737.4; %[km] Earth polar radius (NASA Fact Sheet)

% Get the points on the ellipsoid defined by the two radii:
[Moonx,Moony,Moonz] = ellipsoid(0,0,0,RM_eq,RM_eq,...
    RM_pol,npts);

% Position earth at desired location 
Moonx = Moonx+pos(1);
Moony = Moony+pos(2);
Moonz = Moonz+pos(3);

% Read in the texturemap for the Earth
cdata = imread('Moon_texturemap.jpg');

% Plot the Earth's surface and use the texturemap as the color:
gh = surf(Moonx,Moony,-Moonz,'FaceColor','texturemap','CData',cdata, ...
    'EdgeColor','none');
axis equal