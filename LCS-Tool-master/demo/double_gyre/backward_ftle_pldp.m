%% attracting FLTE modified from double gyre demo oc LCS Tool Master

load('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/CODAR_filled.mat');
% load('ocean_geostrophic_velocity.mat')
vLat = CODAR_filled.V; % here, I made the CODAR data the same size as the demo data to get the grid interpolation to run with current settings, need to fix this later
vLat = permute(vLat, [3 2 1]);
vLon = CODAR_filled.U;
vLon = permute(vLon, [3 2 1]);

lat = CODAR_filled.lat;
lon = CODAR_filled.lon;
dnum = CODAR_filled.dnum;
% convert uits
time = dnum * 24; % when using datenum, convert to hours
vLat = vLat * 60 * 60 * 24 * (1/4349600); % convert from cm/s to dec deg / day;
vLon = vLon * 60 * 60 * 24 * (1/4349600); % convert from cm/s to dec deg / day;
% vLat = vLat * 60 * 60 * (1/4349600); % convert from cm/s to dec deg / day;
% vLon = vLon * 60 * 60 * (1/4349600); % convert from cm/s to dec deg / day;

%% Input parameters
% T_0 = 1170;
% T_1 = T_0+7;
% timespan = [time(T_0) time(T_1)];
timespan = [0, 10];
epsilon = .1;
amplitude = .1;
omega = pi/5;
min_lat = min(lat);
min_lon = min(lon);
max_lat = max(lat);
max_lon = max(lon);
domain = [min_lon, max_lon; min_lat, max_lat];
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];
% lDerivative = @(t,x,~)derivative(t,x,false,epsilon,amplitude,omega);
interpMethod = 'spline';
vLonInterpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLatInterpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);
lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
incompressible = true;
cgStrainOdeSolverOptions = odeset('relTol',1e-5);
lambda = .93:.01:1.07;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
forceEtaComplexNaN = true;
%%
% Lambda lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [.55,.55;.1,.1];
poincareSection(2).endPosition = [1.53,.45;1.95,.05];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 10*deltaX;
stretchLineOdeSolverOptions = odeset('relTol',1e-6);
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
% Elliptic LCSs
tspan = timespan;
[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'odeSolverOptions',lambdaLineOdeSolverOptions);
ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)];
%
hAxes = setup_figure(domain);
hEllipticLcs = plot_elliptic_lcs(hAxes,ellipticLcs);

%% Hyperbolic attracting LCSs
hAxes = setup_figure(domain);
title(hAxes,'Attracting and elliptic LCSs')

% Copy objects from repelling LCS plot
hEllipticLcs = copyobj(hEllipticLcs,hAxes);
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
stretchLine = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
end

hAxes = setup_figure(domain);
title(hAxes,'Attracting and elliptic LCSs')
% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
drawnow

% % Plot attracting LCSs
% hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchLine,'UniformOutput',false);
% hAttractingLcs = [hAttractingLcs{:}];
% set(hAttractingLcs,'color',attractingColor)
% 
% uistack(hEllipticLcs,'top')