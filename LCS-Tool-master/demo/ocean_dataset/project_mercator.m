function project_mercator(ax)
%
%
% USAGE: project_mercator([ax])
%
% Sets the current axes (or axes specified by the optional argument) to an
% approximate mercator projection by setting the 'DataAspectRatio' using the
% following equation:
%
%   [ 1 cos(mean(get(ax,'ylim'))*pi/180) mean(abs(get(ax,'zlim')))*2 ]
%
% ============================================================================
% $RCSfile: project_mercator.m,v $
% $Source: /home/kerfoot/cvsroot/graphics/project_mercator.m,v $
% $Revision: 1.1.1.1 $
% $Date: 2011/06/28 12:29:24 $
% $Author: kerfoot $
% ============================================================================
%

if isequal(nargin,0)
    ax = gca;
end

mean_lat = mean(get(ax, 'YLim'));

set( ax,...
    'DataAspectRatio', [1 cos(mean_lat*pi/180) 1]);
