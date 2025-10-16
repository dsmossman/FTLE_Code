%%load the mask
load('/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/ocean_dataset/PLDP_CODAR_domain.mat');
load('/Volumes/home/jmv208/SWARM_data/SWARM_CODAR.mat');

load('/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/ocean_dataset/palm_grid_mask_80p.mat');
% zeros are good points, MATLAB is start:step:stop
swarm_grid_mask_80p = grid_mask_80p(1:2:end,1:2:end); % size to 101x101
%%
% % this doesn't work, alpha shape function removes points that it thinks
% % are duplicates
% [x,y] = meshgrid(CODAR.lon,CODAR.lat);
% shp = alphaShape(X(:,1),X(:,2));
% tf = inShape(shp,x,y);

%%flag the bad data points
bad_grid=find(swarm_grid_mask_80p==1);

for i=1:length(CODAR.dnum)
    
    u = CODAR.u(:,:,i);
    u(bad_grid) = NaN;
    CODAR.u_domain(:,:,i) = u;
    
    v = CODAR.v(:,:,i);
    v(bad_grid) = NaN;
    CODAR.v_domain(:,:,i) = v;
    
end

%% try replacing all NaNs with zero and feeding into FTLE code

for i =1:length(CODAR.dnum)
    u = CODAR.u_domain(:,:,i);
    u(isnan(u))=0;
    CODAR.u_zeros(:,:,i) = u;
    
    v = CODAR.u(:,:,i);
    v(isnan(v))=0;
    CODAR.v_domain(:,:,i) = v;
    
end
    
