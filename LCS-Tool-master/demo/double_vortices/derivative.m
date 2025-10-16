% derivative Double vortices velocity field
%
% SYNTAX
% derivative_ = derivative(t,position,useEoV,epsilon,amplitude,omega)
%
% INPUT ARGUMENTS
% t: time
% position: [x1;y1;x2;y2;...;xn;yn]
% useEov: logical that controls use of the equation of variation
% epsilon,amplitude,omega: double gyre parameters
%
% REFERENCE
% DOI:10.1016/j.physd.2005.10.007

function derivative_ = derivative(t,x,useEoV,epsilon,gamma)

validateattributes(t,{'double'},{'scalar'})
validateattributes(x,{'double'},{'column'})
% validateattributes(y,{'double'},{'column'})
% Jackie add a y value here ??
% Cannot use validateattributes to check x is a column vector with an even
% number of elements
if mod(numel(x),2)
    error([mfilename,':numelXNotEven'],'numel(x) = %d is not even',numel(x))
end
validateattributes(x,{'double'},{'column'})
% validateattributes(y,{'double'},{'column'})
validateattributes(useEoV,{'logical'},{'scalar'})
validateattributes(epsilon,{'double'},{'scalar'})
validateattributes(gamma,{'double'},{'scalar'})

if useEoV
    idx1 = 1:6:size(x,1)-5;
    idx2 = 2:6:size(x,1)-4;
else
    idx1 = 1:2:numel(x)-1;
    idx2 = 2:2:numel(x);
end
% add in embedded functions here
v = exp(epsilon)/(2*besseli(0,epsilon));
f = @(s) 1-2*v*exp(epsilon*(cos(s)-1));
x_v = 0.5*gamma*exp(-epsilon*(cos(t/gamma))-1)*integral(f,0,t/gamma,'ArrayValued',true);
y_v = exp(epsilon*cos(t/gamma)-1);

derivative_ = nan(size(x));
% add in derivatives here
derivative_(idx1) = -(((idx2-y_v)/((idx2-x_v).^2+(idx2-y_v).^2))-(idx2+y_v)/((idx1-x_v).^2+(idx2+y_v).^2))+v +((epsilon*idx1)/gamma)*sin(t/gamma);
derivative_(idx2) = (idx1-x_v).*((1./((idx1-x_v).^2+(idx2-y_v).^2))-(1./((idx1-x_v).^2+(idx2+y_v).^2)))-((epsilon*idx2)/gamma)*sin(t/gamma);

if useEoV
    % Define terms of the equation of variation
    idx3 = 3:6:size(x,1)-3;
    idx4 = 4:6:size(x,1)-2;
    idx5 = 5:6:size(x,1)-1;
    idx6 = 6:6:size(x,1);
    % add computed partial derivatives here
    dux = ((-2*(y-y_v)*(x-x_v))/((idx1-x_v).^2+(y-y_v).^2).^2 + (2*(y+y_v)*(x-x_v))/((idx1-x_v).^2+(y+y_v).^2).^2) - ((epsilon*idx1)/gamma)*sin(t/gamma);
    duy = ((x-x_v).^2-(idx2-y_v).^2)/((idx2-y_v).^2+(x-x_v).^2).^2 - ((x-x_v).^2-(idx2+y_v).^2)/((idx2+y_v).^2+(x-x_v).^2).^2;
    dvx = idx1*((1/((x-x_v).^2+(y-y_v).^2))-1/((x-x_v).^2+(y+y_v).^2))+(x-x_v)*((-2*(x-x_v))/((idx1-x_v).^2+(y-y_v).^2).^2-(2*(x-x_v))/((idx1-x_v).^2+(y+y_v).^2).^2);
    dvy = (x-x_v)*((-2*(y-y_v))/((x-x_v).^2+(idx2-y_v).^2).^2+(2*(y-y_v))/((x-x_v).^2+(idx2-y_v).^2).^2);
    
    % Perform matrix multiplication manually
    derivative_(idx3) = dux.*x(idx3) + duy.*x(idx5);
    derivative_(idx4) = dux.*x(idx4) + duy.*x(idx6);
    derivative_(idx5) = dvx.*x(idx3) + dvy.*x(idx5);
    derivative_(idx6) = dvx.*x(idx4) + dvy.*x(idx6);
end
