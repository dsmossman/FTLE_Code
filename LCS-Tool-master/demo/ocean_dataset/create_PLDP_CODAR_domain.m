% how I made the domain for the PLDP CODAR
load('/Volumes/home/jmv208/SWARM_data/SWARM_CODAR.mat');
addpath('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/');

filenames = dir('/Volumes/SWARM_CODAR/nc/RU_SWARM_*.nc');

eval(['lon = ncread(''/Volumes/SWARM_CODAR/nc/', filenames(1).name, ''',''lon'');']); %grab lons from codar nc's
eval(['lat = ncread(''/Volumes/SWARM_CODAR/nc/',filenames(1).name, ''',''lat'');']); %grab lons from codar nc's
[X,Y] = meshgrid(lon, lat); %create grid from lons lats
X = X'; %rotate grid
Y = Y'; %rotate grid

start_1 = [datenum(2020, 2, 9, 12, 0, 0)];
len1 = 24;
lw1 = find(CODAR.dnum == start_1);
ind1 = [lw1:(lw1+len1-1)];

all_u1 = CODAR.u(:,:,ind1);
all_v1 = CODAR.v(:,:,ind1);

for j=1:101
    for k=1:101
        
       ind_nan_u=isnan(all_u1(j,k,:));
       ind_nan_v=isnan(all_v1(j,k,:));
       
       if (sum(ind_nan_u)<=12 & (sum(ind_nan_v)<=12))
           u_avg(j,k)=nanmean(all_u1(j,k,:));
           v_avg(j,k)=nanmean(all_v1(j,k,:));
       else
           u_avg(j,k)=NaN;
           v_avg(j,k)=NaN;
           
       end
    end
end

  
mag1 = sqrt(u_avg.^2 + v_avg.^2); % magnitude of current for all files in 'filenames'
title_string = '2020-01-09 12:00:00 GMT--2020-01-13 16:00:00 GMT';
plotAntarctic_CODAR(X,Y,mag1,u_avg,v_avg,title_string) % function will plot and save fig


usercoast_filename = 'COAST_PLDP.mat';
functName = 'ginput';
[X,bI,ptH,lineH] = makeDomainBoundary(usercoast_filename,functName);