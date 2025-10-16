% double check derivation

close all
clear all
addpath('/Users/jveatch/Documents/MATLAB/LCS-Tool-master/');
%% Input parameters
Gamma = 10; % [m^2/s]
eps = 0; % [s^-1]
omega = 0.28; % [s^-1]
d = 1; % [m]

Vv = Gamma/(4*pi*d); 
timespan = [0,1]; % [s]
domain =[-1,1;-2,2]; % [m]
xv = 0.5; % [m]
yv = 1; % [m]

% approximations for creating dimensionless variables
% x/d -> x
% y/d-> y
domain = domain / d;

% xv/d -> xv
% yv/d -> xv
xv = xv/d;
yv = yv/d;

% (eps/omega)->epsilon
epsilon = eps/omega;

% (2*pi*d*Vv)/Gamma->muv
muv = 2*pi*d*Vv/Gamma;

% Gamma/(2*pi*omega*d^2)->gamma
gamma = Gamma/(2*pi*omega*d^2);

% (Gamma*t)/(2*pi*d^2) -> t
timespan = (Gamma*timespan)/(2*pi*d^2);

resolutionX = 50;
[resolutionY,deltaX] = equal_resolution(domain,resolutionX);
resolution = [resolutionX,resolutionY];

%% Velocity definition
useEoV = true;
lDerivative = @(t,x,~)derivative(t,x,useEoV,epsilon,gamma,muv,xv,yv);
incompressible = true;

%% LCS parameters
% Cauchy-Green strain
cgStrainOdeSolverOptions = odeset('relTol',1e-6);

%% Cauchy-Green strain eigenvalues and eigenvectors
[cgEigenvector,cgEigenvalue] = eig_cgStrain(lDerivative,domain,resolution,timespan,'incompressible',incompressible,'odeSolverOptions',cgStrainOdeSolverOptions,'method', 'equationOfVariation', 'coupledIntegration', true);

hAxes = setup_figure(domain);
title(hAxes,'Repelling and elliptic LCSs')
% Plot finite-time Lyapunov exponent
cgEigenvalue2 = reshape(cgEigenvalue(:,2),fliplr(resolution));
ftle_ = ftle(cgEigenvalue2,diff(timespan));
plot_ftle(hAxes,domain,resolution,ftle_);
colormap(hAxes,flipud(gray))
drawnow