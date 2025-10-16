% derivative Double vortices velocity field
%
% SYNTAX
% derivative_ = derivative(t,position,useEoV,epsilon,amplitude,omega)
%
% INPUT ARGUMENTS
% t: time
% position: [x1;y1;x2;y2;...;xn;yn]
% useEov: logical that controls use of the equation of variation
% epsilon, gamma, muv: vortex pair
%
% REFERENCE
% 
% An analytical study of transport, mixing and chaos in an unsteady vortical flow
%  Rom-Kedar, V.; Leonard, A.; Wiggins, S.
% Pub Date: May 1990
% DOI:10.1017/S0022112090000167
%

function derivative_ = derivative(t,x,useEoV,epsilon,gamma,muv,xv,yv)

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
validateattributes(useEoV,{'logical'},{'scalar'})
validateattributes(epsilon,{'double'},{'scalar'})
validateattributes(gamma,{'double'},{'scalar'})
validateattributes(muv,{'double'},{'scalar'})
validateattributes(xv,{'double'},{'scalar'})
validateattributes(yv,{'double'},{'scalar'})

% when assuming that the vortices are not moving, do not use these
% dxvdt = 1/(2*yv)-muv+(epsilon*xv/gamma)*sin(t/gamma);
% dyvdt =  -(epsilon*yv/gamma)*sin(t/gamma);


if useEoV
    idx1 = 1:6:size(x,1)-5;
    idx2 = 2:6:size(x,1)-4;
else
    idx1 = 1:2:numel(x)-1;
    idx2 = 2:2:numel(x);
end

% forcing force
fx =  (epsilon*x(idx1)/gamma).*sin(t/gamma);
fy = -(epsilon*x(idx2) /gamma).*sin(t/gamma);

derivative_ = nan(size(x));

% idx1 and idx2 are the indexes of the x and y respectively. Therefore, to
% denote x as x(idx2) and y as x(idx2)

% dxdt = -[(y-yv)/((x-xv)^2+(y-yv)^2)-(y+yv)/((x-xv)^2+(y+yv)^2)]-muv+(epsilon*x/gamma)*sin(t/gamma);
derivative_(idx1) = -((x(idx2)-yv)./((x(idx1)-xv).^2+(x(idx2)-yv).^2)-(x(idx2)+yv)./((x(idx1)-xv).^2+(x(idx2)+yv).^2))-muv+fx;
% dydt = (x-xv) *[1/((x-xv)^2+(y-yv)^2)-1/((x-xv)^2+(y+yv)^2)]-(epsilon*y/gamma)*sin(t/gamma);
derivative_(idx2) = (x(idx1)-xv).*(1./((x(idx1)-xv).^2+(x(idx2)-yv).^2)-1./((x(idx1)-xv).^2+(x(idx2)+yv).^2))+fy;

if useEoV
    % Define terms of the equation of variation
    idx3 = 3:6:size(x,1)-3;
    idx4 = 4:6:size(x,1)-2;
    idx5 = 5:6:size(x,1)-1;
    idx6 = 6:6:size(x,1);
    
    % forcing force
    dfx =  (epsilon/gamma).*sin(t/gamma);
    dfy = -(epsilon/gamma).*sin(t/gamma);

    
    % TODO
    dux = ((2*(x(idx2)-yv).*(x(idx1)-xv))./(((x(idx2)-yv).^2+(x(idx1)-xv).^2).^2))-((2*(x(idx2)+yv).*(x(idx1)-xv))./(((x(idx2)+yv).^2+(x(idx1)-xv).^2).^2))+dfx;
    duy = (((x(idx1)).^2-(2*x(idx1)*xv)+(xv.^2)-(yv.^2)+(2*x(idx2).*yv)-(yv.^2))./(((x(idx2)-yv).^2+(x(idx1)-xv).^2).^2))-(((x(idx1)).^2-(2*x(idx1).*xv)+(xv.^2)-(yv.^2)+(2*x(idx2).*yv)-(yv.^2))./(((x(idx2)+yv).^2+(x(idx1)-xv).^2).^2));
    dvx = (((-x(idx1)).^2+(2*x(idx1)*xv)-(xv.^2)+(x(idx2).^2)-(2*x(idx2)*yv)+(yv.^2))./(((x(idx2)-yv).^2+(x(idx1)-xv).^2).^2))+((2*(x(idx1)-xv))./(((x(idx2)+yv).^2+(x(idx1)-xv).^2).^2));
    dvy = -((2*(x(idx1)-xv).*(x(idx2)-yv))./(((x(idx2)-yv).^2+(x(idx1)-xv).^2).^2))+((2*(x(idx2+yv)))./(((x(idx2)+yv).^2+(x(idx1)-xv).^2).^2))+dfy;
    
    % Perform matrix multiplication manually
    derivative_(idx3) = dux.*x(idx3) + duy.*x(idx5);
    derivative_(idx4) = dux.*x(idx4) + duy.*x(idx6);
    derivative_(idx5) = dvx.*x(idx3) + dvy.*x(idx5);
    derivative_(idx6) = dvx.*x(idx4) + dvy.*x(idx6);
end
