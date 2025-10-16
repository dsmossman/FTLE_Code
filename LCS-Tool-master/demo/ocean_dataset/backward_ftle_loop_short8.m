%% attracting FLTE modified from double gyre demo of LCS Tool Master
% taking another wack at it 05April2022
% another wack 30June2022
% outfitted for Baffin 23Aug2022

load('/home/jmv208/SWARM_data/CODAR_filled_edges.mat');
CODAR_filled = CODAR_filled_edges;
clear CODAR_filled_edges
vLat = CODAR_filled.V; 
vLat = permute(vLat, [3 2 1]);
vLon = CODAR_filled.U;
vLon = permute(vLon, [3 2 1]);

lat = CODAR_filled.lat;
lon = CODAR_filled.lon;
dnum = CODAR_filled.dnum;
% convert uits
time = dnum * 24; % when using datenum, convert to hours
vLat = vLat * 60 * 60 * 24 * (1/100) * (1/1000) * (1/111.12); % convert from cm/s to dec deg / hour;
vLon = vLon * 60 * 60 * 24 * (1/100) * (1/1000) * (1/(111.12*cos(-64))); % convert from cm/s to dec deg / hour;

addpath('/home/jmv208/LCS-Tool-master/')
%% Input parameters --> edited for PLDP block of code from ocean_dataset demo

min_lat = min(lat);
min_lon = min(lon);
max_lat = max(lat);
max_lon = max(lon);
domain = [min_lon, max_lon; min_lat, max_lat];
resolutionX = 500;
[resolutionY,~] = equal_resolution(domain,resolutionX);
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

% % Stretch lines
% stretchLineMaxLength = 20;
% stretchLineLocalMaxDistance = 10*deltaX;
% stretchLineOdeSolverOptions = odeset('relTol',1e-6);

% % Graphic properties
% repellingColor = 'r';
% attractingColor = 'b';
% ellipticColor = [0,.6,0];

% start_time = datestr(CODAR_filled.dnum(1));
% end_time = datestr(CODAR_filled.dnum(end - 7));
k = 1;
%%
for i=datetime('28-Feb-2020 00:00:00'):hours(1):datetime('7-Mar-2020 00:00:00')
    date_str = datestr(i);
    T_0 = find(CODAR_filled.dnum == datenum(date_str));
    T_1 = T_0+12;
    timespan = [time(T_0), time(T_1)];

    %% Cauchy-Green strain eigenvalues and eigenvectors
    [~,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);


    %% Hyperbolic attracting LCSs

    % FIXME Part of calculations in seed_curves_from_lambda_max are
    % unsuitable/unecessary for stretchlines do not follow ridges of λ₁
    % minima
%     stretchLine = seed_curves_from_lambda_max(stretchLineLocalMaxDistance,stretchLineMaxLength,-cgEigenvalue(:,1),cgEigenvector(:,3:4),domain,resolution,'odeSolverOptions',stretchLineOdeSolverOptions);

%     hAxes = setup_figure(domain);
%     str = strcat('Attracting and elliptic LCSs ', date_str);
%     title(hAxes,str)
    % Plot finite-time Lyapunov exponent
    cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
    ftle_ = ftle(cgEigenvalue2,diff(timespan));
    ftle_all.ftle(:,:,k) = ftle_;
    ftle_all.time(k) = i;
%     plot_ftle(hAxes,domain,resolution,ftle_);
%     colormap(hAxes,flipud(gray))
%     drawnow
%     hold on;

    % % Plot attracting LCSs
    % hAttractingLcs = cellfun(@(position)plot(hAxes,position(:,1),position(:,2)),stretchLine,'UniformOutput',false);
    % hAttractingLcs = [hAttractingLcs{:}];
    % set(hAttractingLcs,'color',attractingColor)
    % 
    % uistack(hEllipticLcs,'top')
    %% plot Antarctic Map

%     addpath '/home/jmv208/antarcticaPlotting/'
%     addpath '/home/jmv208/antarcticaPlotting/functions/'

%             % the shape file used for this code is split into three different
%         % segments, seperated by NaNs. Older versions of MATLAB do not have a
%         % problem with this, newer versions (like MATLAB_R2019) do not like
%         % this. The code below splits the shape file into three readable
%         % pieces.
%         tanLand = [240,230,140]./255;
%         S1 = shaperead('cst00_polygon_wgs84.shp');
%         S2=S1(1:1174);
%         ind=[0,find(isnan(S1(1175).X))];
%         for x=1:length(ind)-1
%             S2(1174+x)=S1(1175);
%             S2(1174+x).X=S2(1174+x).X(ind(x)+1:ind(x+1));
%             S2(1174+x).Y=S2(1174+x).Y(ind(x)+1:ind(x+1));
%         end
%         mapshow(S2,'facecolor', tanLand)
%     %     set(gca,'color',[ .21, .50, .85]);
%         hold on
% 
%         %marking location of CODAR stations
%         s(1) = plot(-64.0554167, -64.7741833, 'g^',...
%             'markersize', 12,...
%             'markerfacecolor', 'green',...
%             'markeredgecolor', 'black');
%         s(2) = plot(-64.3604167, -64.7871667, 'gs',...
%             'markersize', 12,...
%             'markerfacecolor', 'green',...
%             'markeredgecolor', 'black');
%         s(3) = plot(-64.0446333, -64.9183167, 'gd',...
%             'markersize', 12,...
%             'markerfacecolor', 'green',...
%             'markeredgecolor', 'black');
% 
%         ylim([-65.05 -64.7])
%         xlim([-64.7 -63.8])
%         caxis([0.5 0.7]);
%     %     a=narrow(-63.6928384831366,-64.683125,.3); % place north facing arrow on upper right corner of map
%         l = legend([s], 'Palmer Station', 'Joubin Islands', 'Wauwermans Islands','Location', 'SouthEast');
% 
%         project_mercator;
%         str_dir = strcat('ftle_7hrs_plots/',str, '.png');
%         saveas(figure(1), str_dir);
% %         F(k) = getframe(figure(1));
%         close all;
        k = k+1; %counter
        fprintf('\n hour complete');
        
end

% ftle_all.hAxes = hAxes;
ftle_all.domain = domain;
ftle_all.resolution = resolution;

save('ftle_all_feb28mar7_HR_12' , 'ftle_all');

exit
% fig = figure;
% movie(fig,F,1)
% movie(fig,F,1)

