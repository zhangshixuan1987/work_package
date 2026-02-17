function [wind_mean, wind_variance, wind_40_year] = get_wind_matrix( wind_maxtrix_raw )
%This is a function that returns a 
%                       42by 21 by  12 array of wind speed mean
%                       42by 21 by  12 array of wind speed variance
%                       42by 21 by 552 array of wind speed record
%(averaged to each of 12 month)

wind_mean = zeros(41,21,12);
wind_variance = zeros(41,21,12);
nyear = size(wind_maxtrix_raw,3)/12;
wind_40_year = zeros(41,21,12,nyear);
%wind_maxtrix_raw = netcdf.getVar(nc_f, 6,[0 0 0 0], [41 21 1 552]);
for ts = 1:12
    for lon = 1:41
        for lat = 1:21
            temp = zeros(nyear,1);
            for year = 1:nyear
                temp(year,1) = wind_maxtrix_raw(lon, lat, ts+12*(year-1));
                wind_40_year(lon, lat, ts, year)=temp(year,1);
            end
            
            wind_mean(lon, lat, ts) = mean(temp,1);
            wind_variance(lon, lat, ts)=var(temp,1);
            
        end
    end
end

end

