%% attracting FLTE modified from double gyre demo of LCS Tool Master
% taking another wack at it 05April2022
% another wack 30June2022
% I just crasehd my MATLAB, outfitted for baffin 16April2024 & added some memory saving strategies

% load('/Users/jveatch/Documents/MATLAB/MAB/data/test_filled/MARACOOS_summer2015test.mat');
load('/home/jmv208/LCS-Tool-master/demo/ocean_dataset/MARACOOS_summer2015test.mat');

date_str = '26-June_2015 00:00:00';
indDay = find(MARACOOS_filled.dnum >= datenum(date_str) & MARACOOS_filled.dnum < datenum(date_str) +1);

vLat = MARACOOS_filled.V(:,:,indDay); 
vLat = permute(vLat, [3 1 2]);
vLon = MARACOOS_filled.U(:,:,indDay);
vLon = permute(vLon, [3 1 2]);

lon = MARACOOS_filled.lat;
lat = MARACOOS_filled.lon;
dnum = MARACOOS_filled.dnum(indDay);

clear MARACOOS_filled

% convert uits
time = dnum * 24; % when using datenum, convert to hours
vLat = vLat * 60 * 60 * 24 * (1/100) * (1/1000) * (1/111.12); % convert from cm/s to dec deg / hour;
vLon = vLon * 60 * 60 * 24 * (1/100) * (1/1000) * (1/(111.12*cos(-64))); % convert from cm/s to dec deg / hour;

% addpath('/Users/jveatch/Documents/MATLAB/LCS-Tool-master')
addpath('/home/jmv208/LCS-Tool-master/')
%% Input parameters --> edited for PLDP block of code from ocean_dataset demo
T_0 = find(dnum == datenum(date_str));
T_1 = T_0+23;
timespan = [time(T_0), time(T_1)];

% adjust to make square --> gridpoints will come out square
min_lat = min(lat);
min_lon = min(lon);
max_lat =max(lat);
max_lon = max(lon);
%     ylim([-65.05 -64.7])
%     xlim([-64.7 -63.8])
domain = [min_lon, max_lon; min_lat, max_lat];
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];
fprintf('\n input parameters successful');

%% Velocity definition --> replaced with block of code from ocean_dataset demo
% Units of variables in .mat
% lon, lat   : degree
% time       : hour
% vLat, vLon : degree/hour
% load('ocean_geostrophic_velocity.mat')
interpMethod = 'spline';
vLonInterpolant = griddedInterpolant({time,lat,lon},vLon,interpMethod);
vLatInterpolant = griddedInterpolant({time,lat,lon},vLat,interpMethod);
% vLonInterpolant = griddedInterpolant({time,lon,lat},vLon,interpMethod);
% vLatInterpolant = griddedInterpolant({time,lon,lat},vLat,interpMethod);

lDerivative = @(t,x,~)derivative(t,x,vLonInterpolant,vLatInterpolant);
incompressible = true;
fprintf('\n velocity has been defined');
%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);

% % Lambda lines
% poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
% poincareSection(1).endPosition = [.55,.55;.1,.1];
% poincareSection(2).endPosition = [1.53,.45;1.95,.05];
% [poincareSection.numPoints] = deal(100);
% nPoincareSection = numel(poincareSection);
% for i = 1:nPoincareSection
%     rOrbit = hypot(diff(poincareSection(i).endPosition(:,1)),diff(poincareSection(i).endPosition(:,2)));
%     poincareSection(i).orbitMaxLength = 2*(2*pi*rOrbit);
% end
% lambda = .93:.01:1.07;
% lambdaLineOdeSolverOptions = odeset('relTol',1e-6);
% forceEtaComplexNaN = true;
% 
% % Shrink lines
% shrinkLineMaxLength = 20;
% shrinkLineLocalMaxDistance = 2*deltaX;
% shrinkLineOdeSolverOptions = odeset('relTol',1e-6);
% 
% % Stretch lines
% stretchLineMaxLength = 20;
% stretchLineLocalMaxDistance = 10*deltaX;
% stretchLineOdeSolverOptions = odeset('relTol',1e-6);
% 
% % Graphic properties
% repellingColor = 'r';
% attractingColor = 'b';
% ellipticColor = [0,.6,0];

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);


%% Hyperbolic attracting LCSs

% FIXME Part of calculations in seed_curves_from_lambda_max are
% unsuitable/unecessary for stretchlines do not follow ridges of λ₁
% minima
% stretchLine = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

hAxes = setup_figure(domain);
str = strcat('Attracting and elliptic LCSs ', date_str);
title(hAxes,'Attracting and elliptic LCSs')
% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray));
drawnow
hold on;
caxis([0.1, 0.2]);

print(fullfile('/home/jmv208/LCS-Tool-master/demo/ocean_dataset', 'MAB_ftle_test4.png'), '-dpng', '-r300');
ftle_mab.domain = domain;
ftle_mab.resolution = resolution;
ftle_mab.ftle = ftle_;

save('mab_ftle_test4' , 'ftle_mab');

% % Plot attracting LCSs
% hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchLine,'UniformOutput',false);
% hAttractingLcs = [hAttractingLcs{:}];
% set(hAttractingLcs,'color',attractingColor)
% 
% uistack(hEllipticLcs,'top')
%% plot Antarctic Map

% addpath '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting/'
% addpath '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting/functions/'
% 
%         % the shape file used for this code is split into three different
%     % segments, seperated by NaNs. Older versions of MATLAB do not have a
%     % problem with this, newer versions (like MATLAB_R2019) do not like
%     % this. The code below splits the shape file into three readable
%     % pieces.
%     tanLand = [240,230,140]./255;
%     S1 = shaperead('cst00_polygon_wgs84.shp');
%     S2=S1(1:1174);
%     ind=[0,find(isnan(S1(1175).X))];
%     for x=1:length(ind)-1
%         S2(1174+x)=S1(1175);
%         S2(1174+x).X=S2(1174+x).X(ind(x)+1:ind(x+1));
%         S2(1174+x).Y=S2(1174+x).Y(ind(x)+1:ind(x+1));
%     end
%     mapshow(S2,'facecolor', tanLand)
% %     set(gca,'color',[ .21, .50, .85]);
%     hold on
%     
%     %marking location of CODAR stations
%     s(1) = plot(-64.0554167, -64.7741833, 'b^',...
%         'markersize', 12,...
%         'markerfacecolor', 'blue',...
%         'markeredgecolor', 'black');
%     s(2) = plot(-64.3604167, -64.7871667, 'bs',...
%         'markersize', 12,...
%         'markerfacecolor', 'blue',...
%         'markeredgecolor', 'black');
%     s(3) = plot(-64.0446333, -64.9183167, 'bd',...
%         'markersize', 12,...
%         'markerfacecolor', 'blue',...
%         'markeredgecolor', 'black');
%      
%     %%Plot the Survey Grid
% %     adgr=plot(adelie_lon, adelie_lat, 'r-','linewidth',[2]);
% %     gegr=plot(gentoo_lon, gentoo_lat, 'g-','linewidth',[2]);
%     
% %     ylim([-65.0 -64.7])
% %     xlim([-64.4 -63.6])
%     ylim([-65.05 -64.7])
%     xlim([-64.7 -63.8])
%     caxis([0.5 0.7]);
% %     a=narrow(-63.6928384831366,-64.683125,.3); % place north facing arrow on upper right corner of map
%     l = legend([s], 'Palmer Station', 'Joubin Islands', 'Wauwermans Islands','Location', 'SouthEast');
% %     l = legend([adgr,gegr],'Adelie Transect', 'Gentoo Transect','Location', 'SouthEast');
% 
%     %l = legend([s], 'Palmer Station', 'Wauwermans Islands', 'Joubin Islands','Location', 'SouthEast');
%    
%     project_mercator;

exit
