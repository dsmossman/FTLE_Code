%% Convert ftle data to .nc
% Jackie Veatch, 24July2024

% load in ftle and filled hfr data, mask ftle to exclude edges, save as .nc

addpath(genpath('C:/Users/dmossman/Box/FTLE Work/'))

clear variables; close all;

cd 'C:\Users\dmossman\Box\FTLE Work\Processed Data'

% for year = 2013:2022
    % Load data for the current year
    ftle_file = 'Glider Deployment Data\MARACOOS_spring_deployment_2023_filled.mat';
    maracoos_file = 'Gapfilled Data\MARACOOS_spring_deployment_2023_filled.mat';
    
    load(ftle_file)
    load(maracoos_file)

    % Define FTLE lat and lon
    maxLon = ftle_mab.domain(2,2);
    minLon = ftle_mab.domain(2,1);
    resolutionLon = ftle_mab.resolution(2);
    maxLat = ftle_mab.domain(1,2);
    minLat = ftle_mab.domain(1,1);
    resolutionLat = ftle_mab.resolution(1);

    stepLon = (maxLon - minLon) / (resolutionLon - 1);
    lon_ftle = minLon:stepLon:maxLon;
    stepLat = (maxLat - minLat) / (resolutionLat - 1);
    lat_ftle = minLat:stepLat:maxLat;
    [lon_ftleGrid, lat_ftleGrid] = meshgrid(lon_ftle, lat_ftle);

    lon = double(MARACOOS_filled.lon);
    lat = double(MARACOOS_filled.lat);

    dt = datetime(ftle_mab.time, 'ConvertFrom', 'datenum');

    spring_mask = ismember(month(dt), [4, 5, 6]);
    datePartArray = dateshift(dt(spring_mask), 'start', 'day');
    uniqueDates_F = unique(datePartArray);
    numspringDay = numel(uniqueDates_F);

    ftle_mab.time = ftle_mab.time(spring_mask);
    ftle_mab.ftle = ftle_mab.ftle(:,:,spring_mask);

    % Process FTLE data
    for i = 1:length(ftle_mab.time)
        cvrg = MARACOOS_filled.hfrcvrg(:,:,i);
        cvrg = double(cvrg);

        outside = computeMask(lon, lat, cvrg', 1, lon_ftleGrid, lat_ftleGrid);

        ftleDay = ftle_mab.ftle(:,:,i)';
        ftleDay(outside) = NaN;
        
        ftle_mab.ftle_noEdges(:,:,i) = ftleDay;
    end

    ftle_mab.ftle = permute(ftle_mab.ftle, [2,1,3]);

    time_readable = datestr(ftle_mab.time);

    % Create a new NetCDF file
    outputFile = 'Glider Deployment Data\MARACOOS_spring_deployment_2023_filled.nc';
    ncid = netcdf.create(outputFile, 'CLOBBER');

    % Define the dimensions of the NetCDF file
    dimid1 = netcdf.defDim(ncid, 'Time', size(ftle_mab.time,2));
    dimid2 = netcdf.defDim(ncid, 'Lat', size(lat_ftle,2));
    dimid3 = netcdf.defDim(ncid, 'Lon', size(lon_ftle,2));
    dimid4 = netcdf.defDim(ncid, 'Char', 11);

    % Define the variable in the NetCDF file
    varid1 = netcdf.defVar(ncid, 'ftle', 'double', [dimid2, dimid3, dimid1]);
    varid2 = netcdf.defVar(ncid, 'ftle_noEdges', 'double', [dimid2, dimid3, dimid1]);
    varid3 = netcdf.defVar(ncid, 'time', 'double', [dimid1]);
    varid4 = netcdf.defVar(ncid, 'lon', 'double', [dimid3]);
    varid5 = netcdf.defVar(ncid, 'lat', 'double', [dimid2]);
    varid6 = netcdf.defVar(ncid, 'time_readable', 'char', [dimid1, dimid4]);

    % Add metadata as suggested by NetCDF Climate and Forecast (CF) Metadata Conventions
    netcdf.putAtt(ncid, varid1, 'long_name', 'Finite Time Lyapunov Exponent (1/hr)');
    netcdf.putAtt(ncid, varid1, 'standard_name', 'FTLE');
    netcdf.putAtt(ncid, varid1, 'units', '1/hr');
    netcdf.putAtt(ncid, varid1, 'FillValue', 'NaN');

    netcdf.putAtt(ncid, varid2, 'long_name', 'Finite Time Lyapunov Exponent (1/hr) without edges of field');
    netcdf.putAtt(ncid, varid2, 'standard_name', 'FTLE');
    netcdf.putAtt(ncid, varid2, 'units', '1/hr');
    netcdf.putAtt(ncid, varid2, 'FillValue', 'NaN');

    netcdf.putAtt(ncid, varid3, 'long_name', 'Time');
    netcdf.putAtt(ncid, varid3, 'standard_name', 'time');
    netcdf.putAtt(ncid, varid3, 'units', 'days since 2000-01-01 00:00:00');
    netcdf.putAtt(ncid, varid3, 'FillValue', 'NaN');
    netcdf.putAtt(ncid, varid3, 'calendar', 'gregorian');
    netcdf.putAtt(ncid, varid3, 'TimeZone', 'GMT');

    netcdf.putAtt(ncid, varid4, 'long_name', 'Longitude');
    netcdf.putAtt(ncid, varid4, 'standard_name', 'longitude');
    netcdf.putAtt(ncid, varid4, 'units', 'degrees_east');
    netcdf.putAtt(ncid, varid4, 'FillValue', 'NaN');

    netcdf.putAtt(ncid, varid5, 'long_name', 'Latitude');
    netcdf.putAtt(ncid, varid5, 'standard_name', 'latitude');
    netcdf.putAtt(ncid, varid5, 'units', 'degrees_north');
    netcdf.putAtt(ncid, varid5, 'FillValue', 'NaN');

    netcdf.putAtt(ncid, varid6, 'long_name', 'Time');
    netcdf.putAtt(ncid, varid6, 'standard_name', 'time');
    netcdf.putAtt(ncid, varid6, 'units', 'date_string');
    netcdf.putAtt(ncid, varid6, 'FillValue', 'NaN');
    netcdf.putAtt(ncid, varid6, 'calendar', 'gregorian');
    netcdf.putAtt(ncid, varid6, 'TimeZone', 'GMT');

    % End the definition mode
    netcdf.endDef(ncid);

    % Write the data to the NetCDF file
    netcdf.putVar(ncid, varid1, ftle_mab.ftle);
    netcdf.putVar(ncid, varid2, ftle_mab.ftle_noEdges);
    netcdf.putVar(ncid, varid3, ftle_mab.time);
    netcdf.putVar(ncid, varid4, lon_ftle);
    netcdf.putVar(ncid, varid5, lat_ftle);
    netcdf.putVar(ncid, varid6, time_readable);

    % Close the NetCDF file
    netcdf.close(ncid);
% end

function outside = computeMask(lon, lat, cvrg_avg, contourLevel, lon_ftleGrid, lat_ftleGrid)
    % Compute the contour matrix
    C = contourc(lon, lat, cvrg_avg, [contourLevel contourLevel]);

    % Initialize cell arrays to store coordinates
    lonContours = {};
    latContours = {};

    % Parse the contour matrix
    k = 1;
    while k < size(C, 2)
        n = C(2, k);
        contourLevelValue = C(1, k);

        % Check if the contour level matches the desired level
        if contourLevelValue == contourLevel
            lonContours{end+1} = C(1, k+1:k+n);
            latContours{end+1} = C(2, k+1:k+n);
        end

        k = k + n + 1;
    end

    % Initialize the mask for points outside the contour
    outside = true(size(lon_ftleGrid));

    % Loop through each contour segment to update the mask
    for i = 1:length(lonContours)
        [in, on] = inpolygon(lon_ftleGrid, lat_ftleGrid, lonContours{i}, latContours{i});
        insideOrOn = in | on;
        outside = outside & ~insideOrOn; % Update the mask
    end
end
