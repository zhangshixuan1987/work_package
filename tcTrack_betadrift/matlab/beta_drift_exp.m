clear variables; clc; close all
%load wind data
load('wind_data_exp')

  %% calculate the beta drift for exp
%average over june to november, and using a ratio to sum vertically
alfa = 0.8;
large_scale_wind_4x(:,:,1) = alfa * nanmean(u850_mean(:,:,6:11),3) + (1-alfa) * nanmean(u200_mean(:,:,6:11),3);
large_scale_wind_4x(:,:,2) = alfa * nanmean(v850_mean(:,:,6:11),3) + (1-alfa) * nanmean(v200_mean(:,:,6:11),3);

% find the environmental variables that goes into beta drift regression
du_dy = nan(41,21);
dv_dx = nan(41,21);
f =  nan(41,21);
du_dp = nan(41,21);
dv_dp = nan(41,21);
rela_vort_dy = nan(41,21);
lon_for_reg = nan(41,21);
u850_for_reg = nan(41,21);
v850_for_reg = nan(41,21);

% only calculate for the inside grid, keep the boundaries as nan.
for i = 2:40
    for j =2:20
        lat = lat1(j);
        lon = lon1(i);
        degree_per_meter(lat, 'lon');
        degree_per_meter(lat, 'lat');
        du_dy(i,j) = (large_scale_wind_4x(i,j+1,1)-large_scale_wind_4x(i,j-1,1))/...
            5.0/degree_per_meter(lat, 'lat')*10^(-5);  %UNIT 10^5
        dv_dx(i,j) = (large_scale_wind_4x(i+1,j,2)-large_scale_wind_4x(i-1,j,2))/...
            5.0/degree_per_meter(lat, 'lon')*10^(-5);  %UNIT 10^5
        f(i,j) = 2*7.29*10^(-5)*sin(lat/180*pi)*10^5;  %UNIT 10^-5
        du_dp(i,j) = (mean(u850_mean(i,j,6:11),3)-mean(u200_mean(i,j,6:11),3))/650*1000; %UNIT 10^-3
        dv_dp(i,j) = (mean(v850_mean(i,j,6:11),3)-mean(v200_mean(i,j,6:11),3))/650*1000; %UNIT 10^-3
        lon_for_reg(i,j) = lon;
        u850_for_reg(i,j) = mean(u850_mean(i,j,6:11),3);
        v850_for_reg(i,j) = mean(v850_mean(i,j,6:11),3);
    end
end

for i = 2:40
    for j =2:20
        lat = lat1(j);
        rela_vort_dy(i,j) = (dv_dx(i,j+1) -du_dy(i,j+1)-dv_dx(i,j-1)+du_dy(i,j-1))/...
            5.0/degree_per_meter(lat, 'lat')*10^(-4);  %UNIT 10^9
    end
end

beta_drift_4x(:,:,1) =2.567 * dv_dx +...
      0.9285 * f +...
      0.0873 * du_dp +...
      0.2107 * dv_dp +...
     -0.365  * rela_vort_dy +...
     -0.0195 * lon_for_reg +...
     -0.5766 * u850_for_reg +...
     -0.6692 * v850_for_reg +...
     -0.4534;
 
beta_drift_4x(:,:,2) =-0.0998 * du_dp +...
      0.2493 * dv_dp +...
      0.1749 * u850_for_reg +...
     -1.2009 * v850_for_reg +...
      1.2707;
  

%% save files
save('beta_drift_con', 'beta_drift_con')
save('beta_drift_4x', 'beta_drift_4x')
save('ASO_beta_drift_processing')