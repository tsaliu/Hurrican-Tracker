clc; clear; close all;

worldmap('world')
load coastlines
plotm(coastlat,coastlon)

% load coast
axesm('ups')% include other axesm options as needed to define lat lon range

framem on;
axis off;    % Turns off surrounding box
tightmap

latd = 30;
lond = -50;
pcolorm(latd,lond,40) % pseudocolor plot "stretched" to the grid
% geoshow(lat,long,'Color','k') %add continental outlines 
