%%% plot the particles_all structure (location of drifters from hourly
%%% trajectories) at a certain timestamp specified by index k
        % datestr(particle_all(1).TimeStamp)
%%% this code will also calculate the relative particle density of 'bins'
%%% within the CODAR field, and plot this on a pcolor map behind the location
%%% of drifters.

% INPUT: particle_all (created by 'save_particle_positions.m')
% OUTPUT: a figure of relative particle density and drifter locations at
% TimeStamp k

% the purpose of this figure is to check for the completion of particle_all
% structure and calculate/plot relative particle density. takes a few
% minutes to run.


addpath(genpath('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting'))

%%load the particle densities
load('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/particles_all_penguin_extendedGrid.mat')

%%Calculate the 2D histogram
%%Configure the grid for the 2D histogram
Xmax = -63.9;
Xmin = -64.6;
Ymax = -64.75;
Ymin = -65;

%%calculate the resolution of the bins 
lon_dist=(Xmax-Xmin)*111.12*cos(-64*pi/180);
lat_dist=(Ymax-Ymin)*111.12;

%%RES of the box in km
RES=1;
nbins_y=ceil(lat_dist/RES);
nbins_x=ceil(lon_dist/RES);
Xedges = [-Inf linspace(Xmin,Xmax,nbins_x) Inf];
Yedges = [-Inf linspace(Ymin,Ymax,nbins_y) Inf];
count=[];
count_dens=[];


%%%Test with a certain time 
%k=630;  %%Jan 27, 2015  11:00 GMT
% k=142;  %%Jan 27, 2015  08:00 GMT


for k=1:length(particles_all)
    
    X = particles_all(k).Lon;
    Y = particles_all(k).Lat;
    
    [count(:,:,k),Xedges,Yedges] = histcounts2(X,Y,Xedges,Yedges,'Normalization','count');
    [count_dens(:,:,k),Xedges,Yedges] = histcounts2(X,Y,Xedges,Yedges,'Normalization','countdensity');
    
    median_count=median(reshape(count(6:36,2:27,k), 26*31,1));
    RPD.rpd(:,:,k) = count(6:36,2:27,k)-median_count;
    % h1 = histogram2(X,Y,Xedges,Yedges,'Normalization','count');
    % h1.FaceColor = 'flat';
    % h1.DisplayStyle = 'tile';
    % view(2)
    
end

RPD.lon = Xedges(6:36);
RPD.lat = Yedges(2:27);
RPD.TimeStamp = [particles_all.TimeStamp].';
RPD.createdWith = 'plot_particles_all_RPD.m';
RPD.creadedInDir = pwd;
RPD.createdFor = 'extended grid to include more penguin foraging locations';

%% calculate the scale factor to go from square degrees to square km

% scale_deg_y=(Ymax-Ymin)/nbins_y;
% scale_deg_x=(Xmax-Xmin)/nbins_x;
% area_deg=scale_deg_y*scale_deg_x;
% area_km=(lon_dist/nbins_x)*(lat_dist/nbins_y);
% 
% scale_factor=area_km/area_deg;
%%plot(squeeze(count_dens(8,7,:))/scale_factor)


X_Count=[];
Y_Count=[];

for xx=7:1:37  %%Only include grid points in the HF footprint
    X_Count=[X_Count; mean([Xedges(xx) Xedges(xx-1)])];
end
for yy=3:1:25
    Y_Count=[Y_Count; mean([Yedges(yy) Yedges(yy-1)])];
end

[XX, YY]=meshgrid(X_Count, Y_Count);

%%Standardize based on the average count across the domain

median_count=median(reshape(count(6:34,2:24,k), 23*29,1));


% pc=pcolor(XX,YY,(count(6:34,2:24,k)/median_count)');
% caxis([0 20])
% hold on

%% plot antarctic bathymetry
bathy=load ('antarctic_bathy_2.mat');
ind2= bathy.depthi==99999;
bathy.depthi(ind2)=[];
bathylines1=0:-10:-100;
bathylines2=0:-200:-1400;
bathylines=[bathylines2];

[cs, h1] = contour(bathy.loni,bathy.lati, bathy.depthi,bathylines, 'linewidth', 1);
clabel(cs,h1,'fontsize',6);
set(h1,'LineColor',[0.7 0.7 0.7])
hold on

%%Plot the map background

tanLand = [240,230,140]./255;
S1 = shaperead('antarctica_shape/cst00_polygon_wgs84.shp');
mapshow([S1.X], [S1.Y], 'DisplayType', 'polygon', 'facecolor', tanLand)
hold on

s(1) = plot(-64.0554167, -64.7741833, 'g^',...
    'markersize', 12,...
    'markerfacecolor', 'green',...
    'markeredgecolor', 'black');
s(2) = plot(-64.3604167, -64.7871667, 'gs',...
    'markersize', 12,...
    'markerfacecolor', 'green',...
    'markeredgecolor', 'black');
s(3) = plot(-64.0446333, -64.9183167, 'gd',...
    'markersize', 12,...
    'markerfacecolor', 'green',...
    'markeredgecolor', 'black');


% % ylim([-65.05 -64.7])
% % xlim([-64.5 -63.95])
%     a=narrow(-63.692838483136/6,-64.683125,.3); % place north facing arrow on upper right corner of map
%%l = legend(s, 'Palmer Station', 'Wauwermans Islands', 'Joubin Islands', 'Location', 'SouthEast');

ylim([-65.0 -64.75])
xlim([-64.5 -63.9])

project_mercator;
set(gcf, 'paperposition', [0 0 11 8.5]);
hold on

%%Plot the drifter locations
s_drift=scatter(particles_all(k).Lon, particles_all(k).Lat, 1,'filled','k', 'MarkerEdgeColor','k');
        
        
            