function [ degree_per_meter ] = degree_per_meter( lat_input, type )
%more accurate way of calculating degree per meter
%  source:https://knowledge.safe.com/articles/725/calculating-accurate-length-in-meters-for-lat-long.html

% rlat = reference latitude in radians
rlat = lat_input * pi() /180;

if type == 'lat'
    % Meters per degree Latitude: 
    meter_per_degree = 111132.92 - 559.82 * cos(2* rlat) + 1.175*cos(4*rlat);
    degree_per_meter = 1.0 / meter_per_degree;
end 

if type =='lon'
    % Meters per degree Longitude: 
    meter_per_degree = 111412.84 * cos(rlat) - 93.5 * cos(3*rlat);
    degree_per_meter = 1.0 / meter_per_degree;
end

end

