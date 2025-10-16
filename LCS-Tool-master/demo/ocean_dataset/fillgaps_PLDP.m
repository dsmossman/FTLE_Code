%% This code will fillgaps of CODAR data that has already been transformed from radials to totals
% code was adapted from Hugh
% (/home/hroarty/codar/MARACOOS/Radials2Totals/totals_toolbox/bin/CODAR_driver_totals_QC.m)
% who adapted his code from Erick Fredj (published method)

% This code will 

% need to install data imaging toolbox to run the 'smoothn' function.
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

load('/Volumes/home/jmv208/SWARM_data/SWARM_CODAR.mat');
conf.OI.BaseDir = '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/hourly_codar_fillgaps';
conf.OI.FilePrefix = 'filled_';
conf.Totals.FileSuffix = '.mat';
conf.MonthFlag = false;

TUV.lat = CODAR.lat;
TUV.lon = CODAR.lon;

ind_time = find(CODAR.dnum >= datenum('9-Jan-2020 00:00:00')); %first day of three HFR sites
dnum = CODAR.dnum(ind_time);
u = CODAR.u(:,:,ind_time);
v = CODAR.v(:,:,ind_time);

[lon,lat] = meshgrid(CODAR.lon, CODAR.lat);
lon = reshape(lon, [length(CODAR.lon)*length(CODAR.lat),1]);
lat = reshape(lat, [length(CODAR.lon)*length(CODAR.lat),1]);
for i=1:length(CODAR.lon)
    for k = 1:length(CODAR.lat)
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

for i=datetime(2020,1,9,0,0,0):hours(1):datetime(2020,3,9,12,0,0) % start analysis when all three HFR sites are running
                
% load data

% make data look like the format that Hugh uses
% dtime = '14-feb-2020 13:00:00';
n = datenum(i);
file_ind = find(CODAR.dnum == n);
TUV.dnum = n;
TUV.U = CODAR.u(:,:,file_ind);
TUV.V = CODAR.v(:,:,file_ind);

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
addpath('/Volumes/home/hroarty/codar/totals_toolbox/totals/toolbox_erick_fredj')
Vs = smoothn({U,V},'robust');

TUVosn.U(hfrcvrg)=Vs{1};
TUVosn.V(hfrcvrg)=Vs{2};

TUVosn.U(isnan(TUVosn.U))=0;
TUVosn.V(isnan(TUVosn.V))=0;
%%-------------------------------------------------------------------------

% Save results
[tdn,tfn] = datenum_to_directory_filename( conf.OI.BaseDir, TUV.dnum, ...
                                           conf.OI.FilePrefix, ...
                                           conf.Totals.FileSuffix, conf.MonthFlag );
tdn = tdn{1};%

if ~exist( tdn, 'dir' )
  mkdir(tdn);
end
save(fullfile(tdn,tfn{1}),'TUV','TUVosn')
% saves file containing origional U, V, lat, lon, time (TUV)
% as well as smoothed U,V data inside domain (TUVosn)

% save to structure
CODAR_filled.U(:,:,k) = TUVosn.U;
CODAR_filled.V(:,:,k) = TUVosn.V;
CODAR_filled.dnum(k) = TUV.dnum;
CODAR_filled.time(k) = CODAR.time(file_ind);
               
k = k+1; %counter
end

CODAR_filled.lon(:,:) = TUVosn.lon;
CODAR_filled.lat(:,:) = TUV.lat;


%% crop 

% CODAR_filled.U = CODAR_filled.U(45:78, 48:81, :);
% CODAR_filled.V = CODAR_filled.V(45:78, 48:81, :);
% CODAR_filled.lon = CODAR_filled.lon(45:78);
% CODAR_filled.lat = CODAR_filled.lat(48:81);