%% attracting FLTE modified from double gyre demo of LCS Tool Master
% taking another wack at it 05April2022
% another wack 30June2022
% I just crasehd my MATLAB, outfitted for baffin 16April2024 & added some memory saving strategies

% load('/Users/jveatch/Documents/MATLAB/MAB/data/test_filled/MARACOOS_summer2015test.mat');
load('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/MARACOOS_summer2015test.mat');

ind = find(MARACOOS_filled.dnum >= datenum('16-Jul_2015 00:00:00'));

vLat = MARACOOS_filled.V(:,:,ind); 
vLat = permute(vLat, [3 1 2]);
vLon = MARACOOS_filled.U(:,:,ind);
vLon = permute(vLon, [3 1 2]);

lon = MARACOOS_filled.lat;
lat = MARACOOS_filled.lon;
dnum = MARACOOS_filled.dnum(ind);

clear MARACOOS_filled

% convert uits
time = dnum * 24; % when using datenum, convert to hours
vLat = vLat * 60 * 60 * 24 * (1/100) * (1/1000) * (1/111.12); % convert from cm/s to dec deg / hour;
vLon = vLon * 60 * 60 * 24 * (1/100) * (1/1000) * (1/(111.12*cos(-64))); % convert from cm/s to dec deg / hour;

% addpath('/Users/jveatch/Documents/MATLAB/LCS-Tool-master')
addpath('/home/jmv208/LCS-Tool-master/')
%% Input parameters --> edited for PLDP block of code from ocean_dataset demo

% adjust to make square --> gridpoints will come out square
min_lat = min(lat);
min_lon = min(lon);
max_lat =max(lat);
max_lon = max(lon);
%     ylim([-65.05 -64.7])
%     xlim([-64.7 -63.8])
domain = [min_lon, max_lon; min_lat, max_lat];
resolutionX = 300;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];
fprintf('\n input parameters successful');

%% Velocity definition --> replaced with block of code from ocean_dataset demo
% Units of variables in .mat
% lon, lat   : degree
% time       : hour
% vLat, vLon : degree/hour
% load('ocean_geostrophic_velocity.mat')
interpMethod = 'spline';
vLonInterpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLatInterpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);
% vLonInterpolant = griddedInterpolant({time,lon,lat},vLon,interpMethod);
% vLatInterpolant = griddedInterpolant({time,lon,lat},vLat,interpMethod);

lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
incompressible = true;
fprintf('\n velocity has been defined');
%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

%% variable to make plot
% Add the land plot
addpath('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/grdread2/')
addpath('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/m_map/')

% Provide the full path to your .grd file
grd_file = '/home/jmv208/LCS-Tool-master/demo/ocean_dataset/GMRTv4_2_20240326topo.grd';
[X, Y, Z] = grdread2(grd_file);

% Plot just the land
mask = Z > 0;
Z_masked = Z;
Z_masked(~mask) = NaN;  % Set values outside mask to NaN

tan_color = [210, 180, 140] / 255;

k = 1;
%% at the start of each day, calculate 24 integration of FTLE
for i=datetime(min(dnum), 'ConvertFrom', 'datenum'):days(1):datetime(max(dnum), 'ConvertFrom', 'datenum')
    date_str = datestr(i);
    T_0 = find(dnum == datenum(date_str));
    T_1 = T_0+23; %% CHANGE INTEGRATION TIME HERE
    timespan = [time(T_0), time(T_1)];

    %% Cauchy-Green strain eigenvalues and eigenvectors
    [~,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

    cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
    ftle_ = ftle(cgEigenvalue2,diff(timespan));
    ftle_mab.ftle(:,:,k) = ftle_;
    ftle_mab.time(k) = datenum(date_str);
    k = k+1;
    str = strcat('\n day complete', date_str);
    fprintf(str);

    % Set up the figure for the first plot
    hAxes = setup_figure(flip(domain));
    str = strcat('Attracting and elliptic LCSs ', date_str);
    title(hAxes, str);

    % Plot finite-time Lyapunov exponent
    plot_ftle(hAxes, flip(domain), flip(resolution), ftle_');
    colormap(hAxes, flipud(gray));
    caxis(hAxes, [0.05, 0.2]); % Set color axis limits for FTLE plot
    hold on;
    
    contourf(X, Y, Z_masked, 0, 'FaceColor', tan_color, 'LineStyle', 'none'); % specify one contour level
    project_mercator;
    
    str = strcat(num2str(k), date_str, '.png');
    print(fullfile('/home/jmv208/LCS-Tool-master/demo/ocean_dataset', str), '-dpng', '-r300');
    
end

ftle_mab.domain = domain;
ftle_mab.resolution = resolution;

save('mab_ftle_testSeason16Jul' , 'ftle_mab');


exit
