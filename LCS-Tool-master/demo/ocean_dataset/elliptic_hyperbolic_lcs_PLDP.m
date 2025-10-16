% % Adapted from Haller's LCS toolbox demo for ocean_dataset
% % must be used with LCS toolbox and the derivative function within the
% % ocean_dataset demo directory
% 
addpath('/Users/jveatch/Documents/MATLAB/LCS-Tool-master'); % path to LCS toolbox

%% load data
load('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Data/my_data/CODAR_filled.mat');
% load('ocean_geostrophic_velocity.mat')
vLat = CODAR_filled.V; % here, I made the CODAR data the same size as the demo data to get the grid interpolation to run with current settings, need to fix this later
vLat = permute(vLat, [3 2 1]);
vLon = CODAR_filled.U;
vLon = permute(vLon, [3 2 1]);

lat = CODAR_filled.lat;
lon = CODAR_filled.lon;
dnum = CODAR_filled.dnum;
% time = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]';
%% convert uits
time = dnum * 24; % when using datenum, convert to hours
vLat = vLat * 60 * 60 * 24 * (1/4349600); % convert from cm/s to dec deg / day;
vLon = vLon * 60 * 60 * 24 * (1/4349600); % convert from cm/s to dec deg / day;
%% Input parameters
% timespan = [100,130];
% timespan = [0,3];
T_0 = find(CODAR_filled.dnum == datenum('01-Mar_2020 00:00:00'));
% T_0 = 1170;
T_1 = T_0+16;
timespan = [time(T_0), time(T_1)];
min_lat = min(lat);
min_lon = min(lon);
max_lat = max(lat);
max_lon = max(lon);
domain = [min_lon, max_lon; min_lat, max_lat];
% domain = [-64.35, -64.1; -64.95, -64.8];
% resolutionX = 1000;
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];
fprintf('\n input parameters successful');
%% Velocity definition
% Units of variables in .mat
% lon, lat   : degree
% time       : hour
% vLat, vLon : degree/hour
% load('ocean_geostrophic_velocity.mat')
interpMethod = 'spline';
vLonInterpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLatInterpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);
lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
incompressible = true;
fprintf('\n velocity has been defined');
%% LCS parameters
% Cauchy-Green strain
cgEigenvalueFromMainGrid = false;
cgAuxGridRelDelta = .01;

% Lambda lines
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
poincareSection(1).endPosition = [3.3,-32.1;3.7,-31.6];
poincareSection(2).endPosition = [1.3,-30.9;2.0,-31.2];
poincareSection(3).endPosition = [4.9,-29.6;5.7,-29.6];
poincareSection(4).endPosition = [4.9,-31.4;5.3,-31.4];
poincareSection(5).endPosition = [3.0,-29.3;3.5,-29.3];
[poincareSection.numPoints] = deal(100);
nPoincareSection = numel(poincareSection);
for i = 1:nPoincareSection
    rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
    poincareSection(i).orbitMaxLength = 4*(2*pi*rOrbit);
end
lambda = .9:.02:1.1;
lambdaLineOdeSolverOptions = odeset('relTol',1e-6,'initialStep',1e-2);
forceEtaComplexNaN = true;

% Shrink lines
shrinkLineMaxLength = 20;
shrinkLineLocalMaxDistance = 2*deltaX;
shrinkLineOdeSolverOptions = odeset('relTol',1e-6);

% Stretch lines
stretchLineMaxLength = 20;
stretchLineLocalMaxDistance = 4*deltaX;
stretchLineOdeSolverOptions = odeset('relTol',1e-6);

% Graphic properties
repellingColor = 'r';
attractingColor = 'b';
ellipticColor = [0,.6,0];
initialPositionMarkerSize = 2;

hAxes = setup_figure(domain);
str = strcat('Repelling and elliptic LCSs',{' '}, datestr(dnum(T_0)));
title(hAxes,str);
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

fprintf('\n LCS parameters successful');
%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'eigenvalueFromMainGrid',cgEigenvalueFromMainGrid,'auxGridRelDelta',cgAuxGridRelDelta);

fprintf('\n Cauchy-Green Strain Tensor and eigenvalues determined, begin plotting');

% Plot repelling finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
% ind_low = find(ftle_ < 0.8);
% ftle_(ind_low) = 0;
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(parula))
caxis([0 1]);
drawnow

%% plot map and codar location
addpath '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting/'
addpath '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting/functions/'

        % the shape file used for this code is split into three different
    % segments, seperated by NaNs. Older versions of MATLAB do not have a
    % problem with this, newer versions (like MATLAB_R2019) do not like
    % this. The code below splits the shape file into three readable
    % pieces.
    tanLand = [240,230,140]./255;
    S1 = shaperead('cst00_polygon_wgs84.shp');
    S2=S1(1:1174);
    ind=[0,find(isnan(S1(1175).X))];
    for x=1:length(ind)-1
        S2(1174+x)=S1(1175);
        S2(1174+x).X=S2(1174+x).X(ind(x)+1:ind(x+1));
        S2(1174+x).Y=S2(1174+x).Y(ind(x)+1:ind(x+1));
    end
    mapshow(S2,'facecolor', tanLand)
%     set(gca,'color',[ .21, .50, .85]);
    hold on
    
    %marking location of CODAR stations
    s(1) = plot(-64.0554167, -64.7741833, 'g^',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(2) = plot(-64.3604167, -64.7871667, 'gs',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(3) = plot(-64.0446333, -64.9183167, 'gd',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
     
    %%Plot the Survey Grid
%     adgr=plot(adelie_lon, adelie_lat, 'r-','linewidth',[2]);
%     gegr=plot(gentoo_lon, gentoo_lat, 'g-','linewidth',[2]);
    
%     ylim([-65.0 -64.7])
%     xlim([-64.4 -63.6])
    ylim([-65.0 -64.7])
    xlim([-64.6 -63.8])
%     a=narrow(-63.6928384831366,-64.683125,.3); % place north facing arrow on upper right corner of map
%     l = legend([s], 'Palmer Station', 'Joubin Islands', 'Wauwermans Islands','Location', 'SouthEast');
%     l = legend([adgr,gegr],'Adelie Transect', 'Gentoo Transect','Location', 'SouthEast');

    %l = legend([s], 'Palmer Station', 'Wauwermans Islands', 'Joubin Islands','Location', 'SouthEast');
   
    project_mercator;
    addpath '/Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/CODE'
    cmocean('thermal');
    caxis([0, 0.5]);

%% Elliptic LCSs
% Plot Poincare sections
hPoincareSection = arrayfun(@(input)plot(hAxes,input.endPosition(:,1),input.endPosition(:,2)),poincareSection,'UniformOutput',false);
hPoincareSection = [hPoincareSection{:}];
set(hPoincareSection,'color',ellipticColor)
set(hPoincareSection,'LineStyle','--')
set(hPoincareSection,'marker','o')
set(hPoincareSection,'MarkerFaceColor',ellipticColor)
set(hPoincareSection,'MarkerEdgeColor','w')
drawnow

[closedLambdaLinePos,closedLambdaLineNeg] = poincare_closed_orbit_range(domain,resolution,cgEigenvector,cgEigenvalue,lambda,poincareSection,'forceEtaComplexNaN',forceEtaComplexNaN,'odeSolverOptions',lambdaLineOdeSolverOptions);

ellipticLcs = elliptic_lcs(closedLambdaLinePos);
ellipticLcs = [ellipticLcs,elliptic_lcs(closedLambdaLineNeg)];

% Plot elliptic LCSs
hEllipticLcs = plot_elliptic_lcs(hAxes,ellipticLcs);
set(hEllipticLcs,'color',ellipticColor)
set(hEllipticLcs,'linewidth',2)

% Plot closed lambda lines
hClosedLambdaLinePos = plot_closed_orbit(hAxes,closedLambdaLinePos);
hClosedLambdaLineNeg = plot_closed_orbit(hAxes,closedLambdaLineNeg);
hClosedLambdaLine = [hClosedLambdaLinePos,hClosedLambdaLineNeg];
set(hClosedLambdaLine,'color',ellipticColor)
drawnow


%% Hyperbolic repelling LCSs
[shrinkLine,shrinkLineInitialPosition] = seed_curves_from_lambda_max(shrinkLineLocalMaxDistance,shrinkLineMaxLength,cgEigenvalue(:,2),cgEigenvector(:,1:2),domain,resolution,'odeSolverOptions',shrinkLineOdeSolverOptions);

% Remove shrink lines inside elliptic LCSs
for i = 1:nPoincareSection
    shrinkLine = remove_strain_in_elliptic(shrinkLine,ellipticLcs{i});
    idx = inpolygon(shrinkLineInitialPosition(1,:),shrinkLineInitialPosition(2,:),ellipticLcs{i}(:,1),ellipticLcs{i}(:,2));
    shrinkLineInitialPosition = shrinkLineInitialPosition(:,~idx);
end

% Plot hyperbolic repelling LCSs
hRepellingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),shrinkLine,'UniformOutput',false);
hRepellingLcs = [hRepellingLcs{:}];
set(hRepellingLcs,'color',repellingColor)
hShrinkLineInitialPosition = arrayfun(@(idx)plot(hAxes,shrinkLineInitialPosition(1,idx),shrinkLineInitialPosition(2,idx)),1:size(shrinkLineInitialPosition,2),'UniformOutput',false);
hShrinkLineInitialPosition = [hShrinkLineInitialPosition{:}];
set(hShrinkLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
set(hShrinkLineInitialPosition,'marker','o')
set(hShrinkLineInitialPosition,'MarkerEdgeColor','w')
set(hShrinkLineInitialPosition,'MarkerFaceColor',repellingColor)

uistack(hEllipticLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
drawnow

%% Hyperbolic attracting LCSs
hAxes = setup_figure(domain);
title(hAxes,'Attracting and elliptic LCSs')
xlabel(hAxes,'Longitude (\circ)')
ylabel(hAxes,'Latitude (\circ)')

% Plot finite-time Lyapunov exponent
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))

% Copy objects from repelling LCS plot
hPoincareSection = copyobj(hPoincareSection,hAxes);
hEllipticLcs = copyobj(hEllipticLcs,hAxes);
hClosedLambdaLine = copyobj(hClosedLambdaLine,hAxes);
drawnow

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
[stretchLine,stretchLineInitialPosition] = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

% Remove stretch lines inside elliptic LCSs
for i = 1:nPoincareSection
    stretchLine = remove_strain_in_elliptic(stretchLine,ellipticLcs{i});
    idx = inpolygon(stretchLineInitialPosition(1,:),stretchLineInitialPosition(2,:),ellipticLcs{i}(:,1),ellipticLcs{i}(:,2));
    stretchLineInitialPosition = stretchLineInitialPosition(:,~idx);
end

% Plot hyperbolic attracting LCSs
hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchLine,'UniformOutput',false);
hAttractingLcs = [hAttractingLcs{:}];
set(hAttractingLcs,'color',attractingColor)
hStretchLineInitialPosition = arrayfun(@(idx)plot(hAxes,stretchLineInitialPosition(1,idx),stretchLineInitialPosition(2,idx)),1:size(stretchLineInitialPosition,2),'UniformOutput',false);
hStretchLineInitialPosition = [hStretchLineInitialPosition{:}];
set(hStretchLineInitialPosition,'MarkerSize',initialPositionMarkerSize)
set(hStretchLineInitialPosition,'marker','o')
set(hStretchLineInitialPosition,'MarkerEdgeColor','w')
set(hStretchLineInitialPosition,'MarkerFaceColor',attractingColor)

uistack(hEllipticLcs,'top')
uistack(hClosedLambdaLine,'top')
uistack(hPoincareSection,'top')
