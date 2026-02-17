%land_or_ocean.m -- return land/ocean for input points
%Purpose: determine if input points are on land or ocean
%
% Syntax:  [isOcean] = land_or_ocean(lat,lon,coastal_res)
%
% Inputs:
%   lat [deg N [-90,90]] - vector of latitude values
%   lon [deg E (-180,180]] - vector of corresponding longitude values
%   coastal_res [pts/deg] - resolution of coastline (gridpts per deg);
%       NOTE: higher resolution = more computing time (1 = decent coarse
%           res; 10 = decent high res)
%
% Outputs:
%   is_Ocean - 1 = ocean, 0 = land
%
% Example: (see land_or_ocean_example.m)
%   lat = -90:10:90;
%   lon = -180:20:180;
%   coastal_res = 5;
%   [isOcean] = land_or_ocean(lat,lon,coastal_res)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: coast.mat (included in MATLAB)

% Original code: http://www.mathworks.com/matlabcentral/answers/1065-determining
%   -whether-a-point-on-earth-given-latitude-and-longitude-is-on-land-or-ocean

% Author: Dan Chavas (adapted from Brett Shoelson 9 Feb 2011 -- see above)
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: http://www.princeton.edu/~dchavas/
% 27 Jan 2014; Last revision: 27 Jan 2014

%------------- BEGIN CODE --------------


function [isOcean] = land_or_ocean_without_load_coast(lat,lon,coastal_res,make_plot, coast)

switch nargin
    case 2
        coastal_res = 1;
        make_plot = 0;
    case 3
        make_plot = 0;
end

if(sum(lon>180)>0)
    lon(lon>180) = lon(lon>180) - 360;  %adjust if using [0,360) lon values
    sprintf('Adjusting lon values from [0,360) to (-180,180]')
end

%% Load coastal data
% coast = load('coast.mat');
%% Define search region (want as small as possible to minimize computation)
lat_search_min = max([min(lat)-2 -90]); %deg N (-90,90]
lat_search_max = min([max(lat)+2 90]);  %deg N (-90,90]
lon_search_min = max([min(lon)-2 -180]);    %deg W (-180,180]
lon_search_max = min([max(lon)+2 180]); %deg W (-180,180]
%% Define land inside of coast
[Z, R] = vec2mtx(coast.lat, coast.long, ...
    coastal_res, [lat_search_min lat_search_max], [lon_search_min lon_search_max], 'filled');

%% Return land/ocean for each input point
val = ltln2val(Z, R, lat, lon);
isOcean = val == 2;
%isLand = ~isOcean;
%% Plot the points on geographic map
if(make_plot)
    figure; worldmap(Z, R)
    geoshow(Z, R, 'DisplayType', 'texturemap')
    colormap([0 1 0;0 0 0;0 1 0;0 0 1])
    plotm(lat,lon,'ro')
end

end

%------------- END OF CODE --------------