% Calculating Lagragian Time Scale for gridpoints in domain of SWARM CODAR
% using autocorrelation function

subname = 'palmer_hfr_acf.m';

global ftbverbose;
if ftbverbose
    fprintf('\nbegin : %s\n', subname)
end

load('/Volumes/home/jmv208/SWARM_data/SWARM_CODAR.mat');
addpath('/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/ocean_dataset');
addpath('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code');
load('swarm_grid_mask_80p.mat');

Domain = 'Palmer';
% Palmer ROI
latlim = [-65.5 -64.5];
lonlim = [-66.5 -63.0];

% READ MATFILE file
times = CODAR.dnum; % in days
numFrames = numel(times);

palmer_lon_rho = CODAR.lon; % decimal degree
palmer_lat_rho = CODAR.lat; % decimal degree
palmer_u_east_surface = CODAR.u; % in cm/s
palmer_v_north_surface = CODAR.v; % in cm/s
domain = find((swarm_grid_mask_80p == 0));
%%
for k = 1:length(domain)
    for i = 1:numFrames
        frame = palmer_u_east_surface(:,:,i);
        values = frame(domain);
        palmer_u_east_surface_selected(i)=values(k);
    end
    test=(isnan(palmer_u_east_surface_selected));
    if sum(test) > 600 % catches bad gridpoints
        noralizedACF_all_u(:,k) = NaN(size(normalizedACF));
        lags_all_u(:,k) = NaN(size(lags));
        unnormalizedACF_all_u(:,k) = NaN(size(unnormalizedACF));
    else
        [normalizedACF, lags] = autocorr(palmer_u_east_surface_selected,'NumLags',1000);
        unnormalizedACF = normalizedACF*var(palmer_u_east_surface_selected,1);
        normalizedACF_all_u(:,k) = normalizedACF;
        lags_all_u(:,k) = lags;
        unnormalizedACF_all_u(:,k) = unnormalizedACF;
    end
    clear palmer_u_east_surface_selected
end

normalizedACF_u = nanmean(normalizedACF_all_u,2);

for i = 1:length(normalizedACF_u)
    SEM = std(normalizedACF_all_u(i,:), 'omitnan')/sqrt(length(normalizedACF_all_u(i,:)));
    ts = tinv([0.025  0.975],length(normalizedACF_all_u(i,:))-1);
    ci = nanmean(normalizedACF_all_u(i,:),2) + ts.*SEM;
    ci_normalized_acf_all_u(i,:) = ci;
    ci_dist = abs(ci-nanmean(normalizedACF_all_u(i,:)));
    ci_dist_normalized_acf_all_u(i,:) = ci_dist;
end

for k = 1:length(domain)
    for i = 1:numFrames
        frame = palmer_v_north_surface(:,:,i);
        values = frame(domain);
        palmer_v_north_surface_selected(i)=values(k);
    end
    test=(isnan(palmer_v_north_surface_selected));
    if sum(test) > 600 % catches bad gridpoints
        noralizedACF_all_v(:,k) = NaN(size(normalizedACF));
        lags_all_v(:,k) = NaN(size(lags));
        unnormalizedACF_all_v(:,k) = NaN(size(unnormalizedACF));
    else
        [normalizedACF, lags] = autocorr(palmer_v_north_surface_selected,'NumLags',1000);
        unnormalizedACF_v = normalizedACF*var(palmer_v_north_surface_selected,1);
        normalizedACF_all_v(:,k) = normalizedACF;
        lags_all_v(:,k) = lags;
        unnormalizedACF_all_v(:,k) = unnormalizedACF;
    end
    clear palmer_v_north_surface_selected
end

normalizedACF_v = nanmean(normalizedACF_all_v,2);

for i = 1:length(normalizedACF_v)
    SEM = std(normalizedACF_all_v(i,:), 'omitnan')/sqrt(length(normalizedACF_all_v(i,:)));
    ts = tinv([0.025  0.975],length(normalizedACF_all_v(i,:))-1);
    ci = nanmean(normalizedACF_all_v(i,:),2) + ts.*SEM;
    ci_normalized_acf_all_v(i,:) = ci;
    ci_dist = abs(ci-nanmean(normalizedACF_all_v(i,:)));
    ci_dist_normalized_acf_all_v(i,:) = ci_dist;
end

for k = 1:length(domain)
    for i = 1:numFrames
        frame_u = palmer_u_east_surface(:,:,i);
        values_u = frame_u(domain);
        frame_v = palmer_v_north_surface(:,:,i);
        values_v = frame_v(domain);
        mag = sqrt(values_u(k)^2 + values_v(k)^2);
        palmer_mag_surface_selected(i)=mag;
    end
    test=(isnan(palmer_mag_surface_selected));
    if sum(test) > 600 % catches bad gridpoints
        noralizedACF_all_mag(:,k) = NaN(size(normalizedACF));
        lags_all_mag(:,k) = NaN(size(lags));
        unnormalizedACF_all_mag(:,k) = NaN(size(unnormalizedACF));
    else
        [normalizedACF, lags] = autocorr(palmer_mag_surface_selected,'NumLags',1000);
        unnormalizedACF = normalizedACF*var(palmer_mag_surface_selected,1);
        normalizedACF_all_mag(:,k) = normalizedACF;
        lags_all_mag(:,k) = lags;
        unnormalizedACF_all_mag(:,k) = unnormalizedACF;
    end
    clear palmer_mag_surface_selected
end
normalizedACF_mag = nanmean(normalizedACF_all_mag,2);

for i = 1:length(normalizedACF_mag)
    SEM = std(normalizedACF_all_mag(i,:), 'omitnan')/sqrt(length(normalizedACF_all_mag(i,:)));
    ts = tinv([0.025  0.975],length(normalizedACF_all_mag(i,:))-1);
    ci = nanmean(normalizedACF_all_mag(i,:),2) + ts.*SEM;
    ci_normalized_acf_all_mag(i,:) = ci;
    ci_dist = abs(ci-nanmean(normalizedACF_all_mag(i,:)));
    ci_dist_normalized_acf_all_mag(i,:) = ci_dist;
end

%%

for k = 1:length(domain)
    for i = 1:numFrames
        frame_u = palmer_u_east_surface(:,:,i);
        values_u = frame_u(domain);
        frame_v = palmer_v_north_surface(:,:,i);
        values_v = frame_v(domain);
        ang = atan(values_v(k)/values_u(k));
        palmer_ang_surface_selected(i)=ang;
    end
    test=(isnan(palmer_ang_surface_selected));
    if sum(test) > 600 % catches bad gridpoints
        noralizedACF_all_ang(:,k) = NaN(size(normalizedACF));
        lags_all_ang(:,k) = NaN(size(lags));
        unnormalizedACF_all_ang(:,k) = NaN(size(unnormalizedACF));
    else
        [normalizedACF, lags] = autocorr(palmer_ang_surface_selected,'NumLags',1000);
        unnormalizedACF = normalizedACF*var(palmer_ang_surface_selected,1);
        normalizedACF_all_ang(:,k) = normalizedACF;
        lags_all_ang(:,k) = lags;
        unnormalizedACF_all_ang(:,k) = unnormalizedACF;
    end
    clear palmer_mag_surface_selected
end
normalizedACF_ang = nanmean(normalizedACF_all_ang,2);

for i = 1:length(normalizedACF_ang)
    SEM = std(normalizedACF_all_ang(i,:), 'omitnan')/sqrt(length(normalizedACF_all_ang(i,:)));
    ts = tinv([0.025  0.975],length(normalizedACF_all_ang(i,:))-1);
    ci = nanmean(normalizedACF_all_ang(i,:),2) + ts.*SEM;
    ci_normalized_acf_all_ang(i,:) = ci;
    ci_dist = abs(ci-nanmean(normalizedACF_all_ang(i,:)));
    ci_dist_normalized_acf_all_ang(i,:) = ci_dist;
end
%%
errorbar(lags,normalizedACF_u, ci_dist_normalized_acf_all_u(:,1), ci_dist_normalized_acf_all_u(:,2));
hold on
errorbar(lags,normalizedACF_v, ci_dist_normalized_acf_all_v(:,1), ci_dist_normalized_acf_all_v(:,2));
errorbar(lags,normalizedACF_mag, ci_dist_normalized_acf_all_mag(:,1), ci_dist_normalized_acf_all_mag(:,2));
errorbar(lags,normalizedACF_ang, ci_dist_normalized_acf_all_ang(:,1), ci_dist_normalized_acf_all_ang(:,2));
xlabel('Lags in hours');
ylabel('Normalize ACF');
title('averaged autocorrelation function');
yline(0)
legend('east-west','north-south', 'magnitude','angle','zero');
xlim([0,200]);
%%
zero_u = find(normalizedACF_u <0);
LTS_u = lags(zero_u(1));
zero_v = find(normalizedACF_v <0);
LTS_v = lags(zero_v(1));

zero_mag = find(normalizedACF_mag <0);
LTS_mag = lags(zero_mag(1));
zero_ang = find(normalizedACF_ang <0);
LTS_ang = lags(zero_ang(1));

