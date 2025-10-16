% add edges back to filled CODAR to increase residence time
% add data outside of 80p index back to filled data, zero everywhere else
% 27June 2022
% filled data created with 'fillgaps_PLDP', modified from Erick's methods

load('/Volumes/home/jmv208/SWARM_data/SWARM_CODAR.mat');
load '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/CODAR_filled.mat'

% already checked that CODAR_filled is masked the same way at each timestamp
%%
ind_time = find(CODAR.dnum >= datenum('9-Jan-2020 00:00:00')); %first day of three HFR sites
dnum = CODAR.dnum(ind_time);
u = CODAR.u(:,:,ind_time);
v = CODAR.v(:,:,ind_time);


%%
test = CODAR_filled.U(:,:,100);
ind_zero = find(test == 0);

for i = 1:length(CODAR_filled.dnum)
    
    raw = u(:,:,i);
    edges = raw(ind_zero);
    filled = CODAR_filled.U(:,:,i);
    filled(ind_zero) = edges;
    ind_nan = isnan(filled);
    filled(ind_nan) = 0;
    CODAR_filled_edges.U(:,:,i) = filled;

end

for i = 1:length(CODAR_filled.dnum)
    
    raw = v(:,:,i);
    edges = raw(ind_zero);
    filled = CODAR_filled.V(:,:,i);
    filled(ind_zero) = edges;
    ind_nan = isnan(filled);
    filled(ind_nan) = 0;
    CODAR_filled_edges.V(:,:,i) = filled;

end

CODAR_filled_edges.dnum = CODAR_filled.dnum;
CODAR_filled_edges.lon = CODAR_filled.lon;
CODAR_filled_edges.lat = CODAR_filled.lat;


