% This is how I made the coastline file for PLDP in order to create
% boundary of HFR

lon = [-65.1, -64.64];
lat = [-65.1, -64.64];
lon = [-64.6, -63.6];
proj = 'mercator';

fname = '/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/ocean_dataset/COAST_PLDP.mat';
res = 4;
fname = makeCoast(lon,lat,proj,fname,res);