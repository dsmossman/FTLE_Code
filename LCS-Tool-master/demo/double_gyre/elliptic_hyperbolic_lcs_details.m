% Double Vortices
% Current Folder should be "Users/jveatch/Documents/MATLAB/LCS-Tool-master"
addpath(genpath('/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/double_vortices/'));
%% Input parameters
epsilon = .1;
gamma = 0.5;
timespan = [0,10];
domain = [0,2;0,1];
resolutionX = 500;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];
% define cartesian grid for initial conditions of traj, define an auxiliary
% grid for differentiating wrt initial conditions
%% Velocity definition
addpath(genpath('/Users/jveatch/Documents/MATLAB/LCS-Tool-master/demo/double_vortices/'));
lDerivative = @(t,x,~)derivative(t,x,false,epsilon,gamma);
incompressible = true;
% enter dixcrete approximation of flow map into derivative function
% use finite differencing over the auxiliary grip point to compute numerically the derivative of the flow map
%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-5);
% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions);
% compute the CG strain tensorfield, its eigenvalue field, and eigenvector
% fields ov er the initial condition grid

hAxes = setup_figure(domain);
title(hAxes,'Double Vortex FTLE')
% Plot finite-time Lyapunov exponent
% find the maximal stretching
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
% compute the FTLE
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
drawnow

