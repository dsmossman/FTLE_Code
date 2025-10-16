%% attracting FLTE modified from double gyre demo of LCS Tool Master
% taking another wack at it 05April2022
% another wack 30June2022
% I just crasehd my MATLAB, outfitted for baffin 16April2024 & added some memory saving strategies

clear all; close all; clc;

addpath(genpath('C:/Users/Delphine/Box/FTLE Work/'))

cd 'C:/Users/Delphine/Box/FTLE Work/Processed Data/Gapfilled Data/'
fname = uigetfile(".mat","Select the gapfilled data");
load(fname);

% Get specific variables from MARACOOS_filled
lon = double(MARACOOS_filled.lat);
lat = double(MARACOOS_filled.lon);
dnum = MARACOOS_filled.dnum;
u = MARACOOS_filled.U;
v = MARACOOS_filled.V;
hfrcvg = MARACOOS_filled.hfrcvrg;

clear MARACOOS_filled

% convert units
Time = dnum * 24; % when using datenum, convert to hours
% Recall that datenum is days since Jan 0, 0000
% Why? Because Matlab is stupid, that's why

cd 'C:\Users\Delphine\Box\FTLE Work\FTLE_Code\LCS-Tool-master\demo\ocean_dataset'
% Changes working directory to where the necessary functions are
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
fprintf("Input parameters successful\n");

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

%% variable to make plot
% Add the land plot

% Provide the full path to your .grd file
% grd_file = 'C:/Users/Delphine/Box/FTLE Work/FTLE_Code/LCS-Tool-master/demo/ocean_dataset/GMRTv4_2_20240326topo.grd';
% [X, Y, Z] = grdread2(grd_file);
% 
% % Plot just the land
% mask = Z > 0;
% Z_masked = Z;
% Z_masked(~mask) = NaN;  % Set values outside mask to NaN
% 
% tan_color = [210, 180, 140] / 255;
% 
% Delphine note: I don't use this section of the code, so I don't think it is
% necessary to translate into python, but I've left it alone for now
%% Calculate integration of FTLE

uniqueDays = datetime(min(dnum), 'ConvertFrom', 'datenum'):days(1):datetime(max(dnum), 'ConvertFrom', 'datenum');

k = 1;

for i=datetime(min(dnum), 'ConvertFrom', 'datenum'):...
        hours(1):...
        (datetime(max(dnum), 'ConvertFrom', 'datenum')-hours(5)) % for each unit in the timespan of the file (here it's hours)
    
    date_str = datestr(i);
    ind = find(uniqueDays == i);
    
    % if sum(hfrcvg(:,:,ind), 'all') < 3000
    % 
    %     ftle_mab.ftle(:,:,k) = NaN(resolution(2), resolution(1));
    %     ftle_mab.time(k) = datenum(date_str);
    %     k = k+1;
    % 
    %     str = strcat("Day complete, not enough data: ", date_str,"\n");
    %     fprintf(str);
    % 
    % else
    
        T_0 = find(dnum == datenum(date_str));
        T_1 = T_0+5; %% CHANGE INTEGRATION TIME HERE; right now it is 6 hours
        timespan = [Time(T_0), Time(T_1)];
        %%
        indDay = find(dnum >= datenum(date_str) & dnum < datenum(date_str) +1);
        timeDay = Time(indDay);
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
        fprintf("Velocity has been defined\n");

        %% Cauchy-Green strain eigenvalues and eigenvectors
        % vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        [~,cgEigenvalue] = eig_cgStrain_edited(lDerivative,... % <<<<<<<<<<<<<<<<<<<<<<
            domain,resolution,timespan,'incompressible',... % <<<<<<<<<<<<<<<<<<<<<<<<< THE BIG ONE
            incompressible,'odeSolverOptions',cgStrainOdeSolverOptions); % <<<<<<<<<<<<
        % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        sz = size(cgEigenvalue);
        
        % if sz(1) == 101400 % size of complete cgEigenvalue calculation,
        % changes with integration time/size/etc so this isn't a reliable
        % method for whether to say the BC failed or not
            
            cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
            ftle_ = ftle(cgEigenvalue2,diff(timespan));
            ftle_mab.ftle(:,:,k) = ftle_;
            ftle_mab.time(k) = datenum(date_str);
            k = k+1;
            str = strcat("Day complete: ", date_str, "\n");
            fprintf(str);

            % Set up the figure for the first plot
            % hAxes = setup_figure(flip(domain));
            % str = strcat('Attracting and elliptic LCSs ', date_str);
            % title(hAxes, str);

            % Plot finite-time Lyapunov exponent
            % plot_ftle(hAxes, flip(domain), flip(resolution), ftle_');
            % colormap(hAxes, flipud(gray));
            % caxis(hAxes, [0.05, 0.2]); % Set color axis limits for FTLE plot
            % hold on;

            % contourf(X, Y, Z_masked, 0, 'FaceColor', tan_color, 'LineStyle', 'none'); % specify one contour level
            % project_mercator;
            % close;

%             str = strcat(num2str(k), date_str, '.png');
%             print(fullfile('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/MAB_spring15Plots', str), '-dpng', '-r300');

        % else
        %     ftle_mab.ftle(:,:,k) = NaN(resolution(2), resolution(1));
        %     ftle_mab.time(k) = datenum(date_str);
        %     k = k+1;
        % 
        %     str = strcat('\n day complete, BC failed', date_str);
        %     fprintf(str);
        % end
        % 
    % end
end
%%
ftle_mab.domain = domain;
ftle_mab.resolution = resolution;

fname = strcat('C:\Users\Delphine\Box\FTLE Work\Processed Data\Glider Deployment Data\','MARACOOS_',current_season,'_deployment_',yr,'hourly_FTLE.mat');
% fname = strcat('C:\Users\Delphine\Box\FTLE Work\Processed Data\Glider Deployment Data\','MARACOOS_2023-11-13_only_FTLE.mat');
save(fname,"ftle_mab","current_season","yr")
