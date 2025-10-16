%% attracting FLTE modified from double gyre demo of LCS Tool Master
% taking another wack at it 05April2022
% another wack 30June2022
% I just crasehd my MATLAB, outfitted for baffin 16April2024 & added some memory saving strategies

% load('/Users/jveatch/Documents/MATLAB/MAB/data/test_filled/MARACOOS_summer2015test.mat');
% load('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/MARACOOS_summer2015test.mat');
load('/home/jmv208/MAB/MARACOOS_filled_2022_winter.mat');

lon = MARACOOS_filled.lat;
lat = MARACOOS_filled.lon;
dnum = MARACOOS_filled.dnum;
u = MARACOOS_filled.U;
v = MARACOOS_filled.V;
hfrcvg = MARACOOS_filled.hfrcvrg;

clear MARACOOS_filled

% convert uits
time = dnum * 24; % when using datenum, convert to hours

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

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

%% variable to make plot
% % Add the land plot
% addpath('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/grdread2/')
% addpath('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/m_map/')
% 
% % Provide the full path to your .grd file
% grd_file = '/home/jmv208/LCS-Tool-master/demo/ocean_dataset/GMRTv4_2_20240326topo.grd';
% [X, Y, Z] = grdread2(grd_file);
% 
% % Plot just the land
% mask = Z > 0;
% Z_masked = Z;
% Z_masked(~mask) = NaN;  % Set values outside mask to NaN
% 
% tan_color = [210, 180, 140] / 255;

datetimeArray = datetime(dnum, 'ConvertFrom', 'datenum');

% Extract the unique dates
uniqueDays = unique(dateshift(datetimeArray, 'start', 'day'));

k = 1;
%% at the start of each day, calculate 24 integration of FTLE
for i=datetime(min(dnum), 'ConvertFrom', 'datenum'):days(1):datetime(max(dnum), 'ConvertFrom', 'datenum')
    
    date_str = datestr(i);
    ind = find(uniqueDays == i);
    
    if sum(hfrcvg(:,:,ind), 'all') < 3000
        
        ftle_mab.ftle(:,:,k) = NaN(resolution(2), resolution(1));
        ftle_mab.time(k) = datenum(date_str);
        k = k+1;
        
        str = strcat('\n day complete, note enough data', date_str);
        fprintf(str);

    else
    
        T_0 = find(dnum == datenum(date_str));
        T_1 = T_0+23; %% CHANGE INTEGRATION TIME HERE
        timespan = [time(T_0), time(T_1)];
        %%
        indDay = find(dnum >= datenum(date_str) & dnum < datenum(date_str) +1);
        timeDay = time(indDay);
        vLat = v(:,:,indDay); 
        vLat = permute(vLat, [3 1 2]);
        vLon = u(:,:,indDay);
        vLon = permute(vLon, [3 1 2]);

        vLat = vLat * 60 * 60 * 24 * (1/100) * (1/1000) * (1/111.12); % convert from cm/s to dec deg / hour;
        vLon = vLon * 60 * 60 * 24 * (1/100) * (1/1000) * (1/(111.12*cos(-64))); % convert from cm/s to dec deg / hour;

        %% Velocity definition --> replaced with block of code from ocean_dataset demo
        % Units of variables in .mat
        % lon, lat   : degree
        % time       : hour
        % vLat, vLon : degree/hour
        % load('ocean_geostrophic_velocity.mat')
        interpMethod = 'spline';
        vLonInterpolant = griddedInterpolant({timeDay,lat,lon},vLon,interpMethod);
        vLatInterpolant = griddedInterpolant({timeDay,lat,lon},vLat,interpMethod);
        % vLonInterpolant = griddedInterpolant({time,lon,lat},vLon,interpMethod);
        % vLatInterpolant = griddedInterpolant({time,lon,lat},vLat,interpMethod);

        lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
        incompressible = true;
        fprintf('\n velocity has been defined');

        %% Cauchy-Green strain eigenvalues and eigenvectors
        [~,cgEigenvalue] = eig_cgStrain_edited(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);

        sz = size(cgEigenvalue);
        
        if sz(1) == 81000 % size of complete cgEigenvalue calculation
            
            cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
            ftle_ = ftle(cgEigenvalue2,diff(timespan));
            ftle_mab.ftle(:,:,k) = ftle_;
            ftle_mab.time(k) = datenum(date_str);
            k = k+1;
            str = strcat('\n day complete', date_str);
            fprintf(str);

%             % Set up the figure for the first plot
%             hAxes = setup_figure(flip(domain));
%             str = strcat('Attracting and elliptic LCSs ', date_str);
%             title(hAxes, str);
% 
%             % Plot finite-time Lyapunov exponent
%             plot_ftle(hAxes, flip(domain), flip(resolution), ftle_');
%             colormap(hAxes, flipud(gray));
%             caxis(hAxes, [0.05, 0.2]); % Set color axis limits for FTLE plot
%             hold on;
% 
%             contourf(X, Y, Z_masked, 0, 'FaceColor', tan_color, 'LineStyle', 'none'); % specify one contour level
%             project_mercator;
% 
%             str = strcat(num2str(k), '.png');
%             print(fullfile('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/MAB_summer14Plots', str), '-dpng', '-r300');

        else
            ftle_mab.ftle(:,:,k) = NaN(resolution(2), resolution(1));
            ftle_mab.time(k) = datenum(date_str);
            k = k+1;
        
            str = strcat('\n day complete, BC failed', date_str);
            fprintf(str);
        end
        
    end
end

ftle_mab.domain = domain;
ftle_mab.resolution = resolution;

save('mab_ftle_winter2022' , 'ftle_mab');


exit
