%% This code will fillgaps of CODAR data that has already been transformed from radials to totals
% code was adapted from Hugh
% (/home/hroarty/codar/MARACOOS/Radials2Totals/totals_toolbox/bin/CODAR_driver_totals_QC.m)
% who adapted his code from Erick Fredj (published method)
% reformatted to fit CONVERGE grid by Jackie 18July2022

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

load('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/hourly_codar_pldp/Converge_2015.mat');
conf.OI.BaseDir = '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/hourly_codar_fillgaps';
conf.OI.FilePrefix = 'filled_';
conf.Totals.FileSuffix = '.mat';
conf.MonthFlag = false;

TUV.lat = Converge_2015.lat;
TUV.lon = Converge_2015.lon;

dnum = Converge_2015.dnum;
u = Converge_2015.u;
v = Converge_2015.v;

[lon,lat] = meshgrid(Converge_2015.lon, Converge_2015.lat);
lon = reshape(lon, [length(Converge_2015.lon)*length(Converge_2015.lat),1]);
lat = reshape(lat, [length(Converge_2015.lon)*length(Converge_2015.lat),1]);
for i=1:length(Converge_2015.lon)
    for k = 1:length(Converge_2015.lat)
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
for i=datetime(2014,11,16,22,0,0):hours(1):datetime(2015,5,26,0,0,0) % start analysis when all three HFR sites are running
                
% load data

% make data look like the format that Hugh uses
% dtime = '14-feb-2020 13:00:00';
n = datenum(i);
file_ind = find(Converge_2015.dnum == n);
TUV.dnum = n;
TUV.U = Converge_2015.u(file_ind,:,:);
TUV.V = Converge_2015.v(file_ind, :,:);

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
Converge15_filled.U(:,:,k) = squeeze(TUVosn.U);
Converge15_filled.V(:,:,k) = squeeze(TUVosn.V);
Converge15_filled.dnum(k) = TUV.dnum;
Converge15_filled.time(k) = Converge_2015.dnum(file_ind);
               
k = k+1; %counter
end

Converge15_filled.lon(:,:) = TUVosn.lon;
Converge15_filled.lat(:,:) = TUV.lat;


%% crop 

% CODAR_filled.U = CODAR_filled.U(45:78, 48:81, :);
% CODAR_filled.V = CODAR_filled.V(45:78, 48:81, :);
% CODAR_filled.lon = CODAR_filled.lon(45:78);
% CODAR_filled.lat = CODAR_filled.lat(48:81);