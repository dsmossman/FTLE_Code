%% fill MAB data one year at a time

clear variables; close all; clc;
addpath(genpath('C:/Users/dmossman/Box/FTLE Work/'))

cd 'C:/Users/dmossman/Box/FTLE Work/Raw HFR Data/'
[ncFile, ncLocation] = uigetfile("*.nc","Select the HFR data");
% Get information about the NetCDF file
ncInfo = ncinfo([ncLocation, ncFile]);

% Loop through each variable in the NetCDF file
for iVar = 1:length(ncInfo.Variables)
    % Read the data for the current variable
    varName = ncInfo.Variables(iVar).Name;
    data = ncread(ncFile, varName);
    
    % Store the data in a structure
    MARACOOS.(varName) = data;
end


%% This code will fillgaps of CODAR data that has already been transformed from radials to totals
% code was adapted from Hugh
% (/home/hroarty/codar/MARACOOS/Radials2Totals/totals_toolbox/bin/CODAR_driver_totals_QC.m)
% who adapted his code from Erick Fredj (published method)
% re-vamped by Jackie 14July2022 to run on MARACOOS HFR, edited 10April2024

% This code will need to install data imaging toolbox to run the 'smoothn' function.
% MATLAB will prompt you to download it with the error it produces if you
% do not already have it installed

%% This section of code was a contribution from Erick Fredj to gap fill the data
% Masking is a destructive process, any totals current inside the mask will be
% set to zero.

% things I need for this to run:
    % conf.OSN.BestCoverageFile --> mask(n,2); made by domain function
    % TUV.LonLat(n,2) --> first column lon, second column lat
    % TUV.U (totals of u)
    % TUV.V (totals of v)
    % conf.OI.BaseDir, dtime, conf.OI.FilePrefix, conf.Totals.FileSuffix, conf.MonthFlag 
    
% some thoughts: The conf strucutre is a great way to keep track of lots of
                % CODAR data... I haven't yet put in the time for it to work well on my
                % computer, but if future Jackie wants to use this method, check out
                % "CODAR_configuration.m"

                
                
% conf.OI.BaseDir = '/Users/jmv208/MATLAB/MAB/data/filled';
% conf.OI.FilePrefix = 'filled_';
% conf.Totals.FileSuffix = '.mat';
% conf.MonthFlag = false;

TUV.lat = MARACOOS.latitude;
TUV.lon = MARACOOS.longitude;

MARACOOS.dt = datetime(MARACOOS.time, 'ConvertFrom', 'posixtime');

MARACOOS.dnum = datenum(MARACOOS.dt);
dnum = MARACOOS.dnum;
u = MARACOOS.u;
v = MARACOOS.v;

[lon,lat] = meshgrid(MARACOOS.longitude, MARACOOS.latitude);
lon = reshape(lon, [length(MARACOOS.longitude)*length(MARACOOS.latitude),1]);
lat = reshape(lat, [length(MARACOOS.longitude)*length(MARACOOS.latitude),1]);

%% For each day, calculate percent coverage and fill gaps

spring_mask = ismember(month(MARACOOS.dt), [4, 5, 6]);
datePartArray = dateshift(MARACOOS.dt(spring_mask), 'start', 'day');
uniqueDates_Sp = unique(datePartArray);
numSpringDays = numel(uniqueDates_Sp);

summer_mask = ismember(month(MARACOOS.dt), [6, 7, 8]);
datePartArray = dateshift(MARACOOS.dt(summer_mask), 'start', 'day');
uniqueDates_Su = unique(datePartArray);
numSummerDays = numel(uniqueDates_Su);

fall_mask = ismember(month(MARACOOS.dt), [10, 11, 12]);
datePartArray = dateshift(MARACOOS.dt(fall_mask), 'start', 'day');
uniqueDates_F = unique(datePartArray);
numFallDays = numel(uniqueDates_F);

winter_mask = ismember(month(MARACOOS.dt), [1, 2, 3]);
datePartArray = dateshift(MARACOOS.dt(winter_mask), 'start', 'day');
uniqueDates_W = unique(datePartArray);
numWinterDays = numel(uniqueDates_W);

seasons = ["spring","summer","fall","winter"];
[numDays, idx] = max([numSpringDays,numSummerDays,numFallDays,numWinterDays]);

current_season = seasons(idx);

uniqueDates = [uniqueDates_Sp,uniqueDates_Su,uniqueDates_F,uniqueDates_W];

c = 1; % counter for saving data to structure
%%
for d = 1:numDays
    % index in time for one day
    indDay = find(MARACOOS.dt >= uniqueDates(d) & MARACOOS.dt < uniqueDates(d)+days(1));
    u = MARACOOS.u(:,:,indDay);
    v = MARACOOS.v(:,:,indDay);
    dnum = MARACOOS.dnum(indDay);

    for i=1:length(MARACOOS.longitude)
        for k = 1:length(MARACOOS.latitude)
            per_cov(i,k)=sum(~isnan(u(i,k,:)))/length(dnum);
                if per_cov(i,k) > 0.8
                    hfrcvrg(i,k) = 1;
                else
                    hfrcvrg(i,k) = 0;
                end
        end
    end
    hfrcvrg = logical(hfrcvrg);
    MARACOOS_filled.hfrcvrg(:,:,d) = hfrcvrg;
%%

    for i=uniqueDates(d):hours(1):uniqueDates(d)+hours(23)
                
        n = datenum(i);
        file_ind = find(MARACOOS.dnum == n);
        if isempty(file_ind)
            MARACOOS_filled.U(:,:,c) = nan(size(TUVosn.U));
            MARACOOS_filled.V(:,:,c) = nan(size(TUVosn.V));
            MARACOOS_filled.dnum(c) = datenum(i);
            fprintf(1, '*********** SKIPPED DAY **********. \n');
        else
            TUV.dnum = n;
            TUV.U = MARACOOS.u(:,:,file_ind);
            TUV.V = MARACOOS.v(:,:,file_ind);

            fprintf(1, 'Starting Smooth Total Field. \n');
            fprintf(1, '---------------------------------------------------------------- \n');

            %% test that inpolygon worked
            % figure
            % 
            % plot(mask(:,1),mask(:,2)) % polygon
            % axis equal
            % hold on
            % 
            % plot(x(hfrcvrg),y(hfrcvrg),'r+') % points inside
            % plot(x(~hfrcvrg),y(~hfrcvrg),'bo') % points outside
            % hold off
            %%
            % Robust Smooth
            TUVosn = TUV;
            TUVosn.CreationInfo= 'Erick Fredj';

            U=TUVosn.U(hfrcvrg);
            V=TUVosn.V(hfrcvrg);

            % set to reset TUVs.U to NaN
            TUVosn.U = NaN(size(TUVosn.U));
            % set to reset TUVs.V to NaN
            TUVosn.V = NaN(size(TUVosn.V));

            %% this function smoothn is located in toolbox_eric_fredj
            Vs = smoothn({U,V},'robust');

            TUVosn.U(hfrcvrg)=Vs{1};
            TUVosn.V(hfrcvrg)=Vs{2};

            TUVosn.U(isnan(TUVosn.U))=0;
            TUVosn.V(isnan(TUVosn.V))=0;
            %%-------------------------------------------------------------------------

    %         % Save results
    %         [tdn,tfn] = datenum_to_directory_filename( conf.OI.BaseDir, TUV.dnum, ...
    %                                                    conf.OI.FilePrefix, ...
    %                                                    conf.Totals.FileSuffix, conf.MonthFlag );
    %         tdn = tdn{1};%
    % 
    %         if ~exist( tdn, 'dir' )
    %           mkdir(tdn);
    %         end
    %         save(fullfile(tdn,tfn{1}),'TUV','TUVosn')
    %         % saves file containing origional U, V, lat, lon, time (TUV)
    %         % as well as smoothed U,V data inside domain (TUVosn)

            % save to structure
            MARACOOS_filled.U(:,:,c) = TUVosn.U;
            MARACOOS_filled.V(:,:,c) = TUVosn.V;
            MARACOOS_filled.dnum(c) = TUV.dnum;
        end
        c = c+1; %counter
    end
    fprintf(1, '*********** DAY FINISHED **********. \n');
    pause(0.5);
    clc;
end
%%
MARACOOS_filled.lon = MARACOOS.longitude;
MARACOOS_filled.lat = MARACOOS.latitude;

yr = string(year(uniqueDates(d)));

fname = strcat('C:\Users\dmossman\Box\FTLE Work\Processed Data\Gapfilled Data\','MARACOOS_',current_season,'_deployment_',yr,'_filled.mat');

save(fname,"MARACOOS_filled","current_season","yr")

% exit