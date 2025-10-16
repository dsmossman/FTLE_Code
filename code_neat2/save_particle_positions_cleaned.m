%%% creates a structure 'particles_all' that contains the position of each
%%% drifter at houlry TimeStamps. The size of the grid on which particles
%%% are released are specified in ' create_hourly_traj_PLDP.m '. The lat-lon grid in this code must
%%% match the particle release grid from the hourly trajectories

%INPUT: hourly trajectory files from ' create_hourly_traj_3Day_PLDP.m '
%OUTPUT: 2D time series of particle positions

% what is happening:
    %searches structure of all Traj files for the same timestamp and puts
    %together a file of the the particle locations for each timestamp
    
% 30Nov2022 cleaned to make a bit more flexible, start_time and end_time
% are still hardcoded

addpath(genpath('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting'))
% addpath(genpath('/Users/kohut/Projects/CONVERGE/PaperProcessing')) % what do I need from here? idk


%%trajectory file directory
traj_dir='/Users/jveatch/Documents/MATLAB/Particle_Track_Code/RPD/code_neat2/penguin_grid/'; %hourly trajectory files

%%Trajectory file directory
filenames = dir([traj_dir,'traj_*.mat']);

start_time=datenum(2020,1,9,0,0,0);
end_time=datenum(2020,2,1,0,0,0);

times=start_time:1/24:end_time;

%integration time cutoff (the first time after release to be used in the 2D
%histogram)
tm_threshold=6; %%Let dirfters move at least six hours from their release time --> gets rid of the "memory" of the grid
tm_cutoff=48; %%Max time to track drifters 
tm_length=(tm_cutoff-tm_threshold)+1;

load([traj_dir, filenames(1).name])

% set vars to correct size for amount of particles released and timesteps
% loaded (TRAJ files loaded)
lon=NaN*ones(length(TRAJ.TrajectoryDuration),tm_length,length(filenames));
latat=NaN*ones(length(TRAJ.TrajectoryDuration),tm_length,length(filenames));
part_times=NaN*ones(tm_length,length(filenames));

%%Setup a loop to integrate drifters from 
for t=1:length(filenames)
    load([traj_dir, filenames(t).name])
    lon(:,:,t)=TRAJ.Lon(:,tm_threshold:tm_cutoff);
    lat(:,:,t)=TRAJ.Lat(:,tm_threshold:tm_cutoff);
    part_times(:,t)=TRAJ.TimeStamp(tm_threshold:tm_cutoff)';
end


%%Build the timeseries of particle positions
for f=1:length(times)
    particles_all(f).Lon=NaN;
    particles_all(f).Lat=NaN;
    particles_all(f).TimeStamp=times(f);
    
    for w=1:length(part_times(1,:))  %%for each release
        
        %ind=find((TRAJ_times(:,w)<=times(f)+(1/72))&(TRAJ_times(:,w)>=times(f)-(1/72)));
        ind=find(part_times(:,w)==times(f));
        
        if ~isempty(ind)
            ind_non_nan=find(~isnan(lon(:,ind,w)));
           
            particles_all(f).Lon=[particles_all(f).Lon; lon(ind_non_nan,ind,w)];
            particles_all(f).Lat=[particles_all(f).Lat; lat(ind_non_nan,ind,w)];
        end
        
    end
end


save /Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/particles_all_penguin_extendedGrid.mat particles_all







