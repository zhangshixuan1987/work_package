function [ A_var, A_cov ] = wind_covariance( x_40_year, y_40_year )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nyear = size(x_40_year,4);
for ts = 1:12
    for lon = 1:41
        for lat = 1:21
            %reshape the data from 4D to [40,1] 40 row 1 column array
            cov_part_1 = reshape(x_40_year(lon, lat, ts,:), [nyear,1]);
            cov_part_2 = reshape(y_40_year(lon, lat, ts,:), [nyear,1]);
            %calculate the covariance. The result is a 2 by 2 array
            cov_matrix = cov(cov_part_1, cov_part_2)
            %the lower triangle of B that satisfy B*B' = cov 
            B = chol(cov_matrix, 'lower');
            A_var (lon, lat, ts) = B(1,1);  
            A_cov (lon, lat, ts) = B(2,1);          
        end
    end
end

end

