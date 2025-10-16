%% This code will fillgaps of CODAR data that has already been transformed from radials to totals
% code was adapted from Hugh
% (/home/hroarty/codar/MARACOOS/Radials2Totals/totals_toolbox/bin/CODAR_driver_totals_QC.m)
% who adapted his code from Erick Fredj (published method)
% re-vamped by Jackie 14July2022 to run on MARACOOS HFR

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
                
addpath('/Users/jveatch/Documents/MATLAB/MAB/m_map'); % newest version of mapping toolbox
addpath('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code');
% load data
root = '/Volumes';
addpath ([ root '/home/codaradm/cocmp_scripts/helpers/']); %location of catTotalStructs
conf.Totals.DomainName='MARA';
% conf.Totals.BaseDir=[root '/home/codaradm/data/totals_maracoos/oi/'];
conf.Totals.BaseDir=[root '/home/codaradm/data/totals/maracoos/oi/mat/5MHz'];
conf.Totals.FilePrefix=['tuv_oi_' conf.Totals.DomainName '_'];
conf.Totals.FileSuffix='.mat';
conf.Totals.MonthFlag=1;
conf.Plot.coastFile = [ '/Users/jveatch/Documents/MATLAB/MAB/njCoast_13Mhz_map_grid.mat'];


starttime = datenum(2022,1,0,0,0,0);
finaltime = datenum(2022,3,31,0,0,0);
tspan = [starttime,finaltime];

 dtime = starttime:1/24:finaltime;
 [f]=datenum_to_directory_filename(conf.Totals.BaseDir,dtime,conf.Totals.FilePrefix,conf.Totals.FileSuffix,conf.Totals.MonthFlag);

%% concatenate the total data
[TUVcat,goodCount] = catTotalStructs(f,'TUV');

ii = false(size(TUVcat)); % I don't really know what this does

%% Put things on a grid
[TUVgrid,gridDim] = gridTotals( TUVcat );

s = size(TUVcat.U);
TUVbest = TUVgrid;
u = reshape( TUVbest.U, [ gridDim, s(2) ] );
v = reshape( TUVbest.V, [ gridDim, s(2) ] );
Lon = reshape( TUVbest.LonLat(:,1), gridDim );
Lat = reshape( TUVbest.LonLat(:,2), gridDim );
dnum = TUVbest.TimeStamp;
%% 
for i=1:185
    for k = 1:155
        per_cov(i,k)=sum(~isnan(u(i,k,:)))/length(dnum);
            if per_cov(i,k) > 0.8
                hfrcvrg(i,k) = 1;
            else
                hfrcvrg(i,k) = 0;
            end
    end
end
hfrcvrg = logical(hfrcvrg);

k = 1;

%%
for i=starttime:1/24:finaltime
                
TUV.dnum = i;
TUV.U = u(:,:,k);
TUV.V = v(:,:,k);

fprintf(1, 'Starting Smooth Total Field. \n');
fprintf(1, '---------------------------------------------------------------- \n');

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
addpath('/Volumes/home/hroarty/codar/totals_toolbox/totals/toolbox_erick_fredj')
Vs = smoothn({U,V},'robust');

TUVosn.U(hfrcvrg)=Vs{1};
TUVosn.V(hfrcvrg)=Vs{2};

TUVosn.U(isnan(TUVosn.U))=0;
TUVosn.V(isnan(TUVosn.V))=0;
%%-------------------------------------------------------------------------

% Save results
[tdn,tfn] = datenum_to_directory_filename( conf.Totals.BaseDir, TUV.dnum, ...
                                           conf.Totals.FilePrefix, ...
                                           conf.Totals.FileSuffix, conf.Totals.MonthFlag );
tdn = tdn{1};%

if ~exist( tdn, 'dir' )
  mkdir(tdn);
end
save(fullfile(tdn,tfn{1}),'TUV','TUVosn')
% saves file containing origional U, V, lat, lon, time (TUV)
% as well as smoothed U,V data inside domain (TUVosn)

% save to structure
MARACOOS_filled.U(:,:,k) = TUVosn.U;
MARACOOS_filled.V(:,:,k) = TUVosn.V;
MARACOOS_filled.dnum(k) = TUV.dnum;
% MARACOOS_filled.time(k) = ;
               
k = k+1; %counter
end

MARACOOS_filled.lat_gridded(:,:) = Lat;
MARACOOS_filled.lon_gridded(:,:) = Lon;


%% crop 

% CODAR_filled.U = CODAR_filled.U(45:78, 48:81, :);
% CODAR_filled.V = CODAR_filled.V(45:78, 48:81, :);
% CODAR_filled.lon = CODAR_filled.lon(45:78);
% CODAR_filled.lat = CODAR_filled.lat(48:81);