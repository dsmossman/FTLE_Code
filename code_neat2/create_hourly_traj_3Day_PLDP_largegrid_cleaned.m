%% Program to release drifters on a fixed grid every hour and tracks them for 3 days
    % Input: is hourly CODAR data from Palmer Deep
    % Output: hourly trajectory files each containing 73 hours of drift data of the particle relesased at the top of that hour
    % to calculate and plot these results, use 'save_particle_positions.m' and then 'plot_particles_all_RPD.m'
    % Jackie adapted this code from Hugh's particle tracks and Erick's TRAJ file creation code
    
    % edits 13June2022 to use larger grid
    % edits 11Nov2022 for flexibility and readability
        % start and stop time in loop are still hardcoded
load '/Users/jveatch/Documents/MATLAB/SWARM/Chapter2/drifter_penguinGrid.mat' % extended grid
% load '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/particle_grid_PLDP_jv_smaller.mat'
% load /Volumes/home/jmv208/SWARM_data/CODAR_filled_edges.mat
load /Users/jveatch/Documents/MATLAB/for_erick/CODAR_filled_edges.mat
% drifter = domain_particle;
    %% generate the drifter release points
    %[wp1]=release_point_generation_NYH; % 121 point grid
%     [wp2]=release_point_generation_matrix_PLDP; % 2500 point grid
    [wp3] = drifter; % 500 point grid (this is the only one that works... seems to be enough data)

    % delete following 4 lines to use smaller grid
    test = wp3 + 0.0075;
    testlong = vertcat(wp3, test);
    lon_tracer = testlong(:,1);
    lat_tracer = testlong(:,2);
    
%     lon_tracer = wp3(:,1);
%     lat_tracer = wp3(:,2);

    ln = lon_tracer;
    lt = lat_tracer;

    LL = [ ln(:), lt(:) ];
    
    % Get some good option values for tracking - otherwise output is junk.
    abs_tol = 1.0e-3; % Not sure about this
    rel_tol = 1.0e-3; % Not sure about this
    maxstep = 1/24; % 1/4 hour  
    options = odeset('RelTol',rel_tol,'AbsTol',abs_tol,'MaxStep',maxstep);
    
% loops through range of dates and releases particles at every hour and
% tracks for three days. Resulting files will contain 73 hours of
% trajectories from each hourly release.

for i=datetime(2020,1,08,12,0,0):hours(1):datetime(2020,2,0,0,0,0)
    
%     DRIFTER=[];
    i = datenum(i);

    dtime=i:1/24:i+3;
    ind = find(CODAR_filled_edges.dnum >= dtime(1) & CODAR_filled_edges.dnum <= dtime(end));
    U = CODAR_filled_edges.U(:,:,ind);
    V = CODAR_filled_edges.V(:,:,ind);
    [Lon, Lat] = meshgrid(CODAR_filled_edges.lon, CODAR_filled_edges.lat);
    
    [year, mon, day, hr, mn, sec]=datevec(i);

    tspan=dtime;
    %drifter=[p.Animation.wp(:,2) p.Animation.wp(:,1)];

    % This calculates trajectories and puts in a TRAJ structure with some
    % extra metadata (but not much).
    TRAJ = ptrack2TRAJstruct('particle_track_ode_grid_LonLat',Lon,Lat,U,V, ...
                         dtime, tspan, ...
                         LL, options );
    TRAJ.TrajectoryDomain = 'PLDP';
    TRAJ.OtherMetadata.ptrack2TRAJstruct.options = options;
    
    filename = strcat('traj_penguinGrid_',num2str(year),'_',num2str(mon),'_',num2str(day),'_',num2str(hr), '_3day.mat');
    filename = fullfile(pwd,filename);
    save(filename,'TRAJ')

end

    