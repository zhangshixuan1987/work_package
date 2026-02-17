% plot

clear all;
clc;
close all;

%% 1. load gensis
load ('sample');

%% 2. load synthetic time series of wind
load('wind_data')

%% plot wind time series
step = 1;
lat = 30;
lon = 270;
time = 270;
plot_wind = nan(48*15-1,15);
while step < 48*15
    u200_mean_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, u200_mean, lat,lon, time);
    v200_mean_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, v200_mean, lat,lon, time);
    u850_mean_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, u850_mean, lat,lon, time);
    v850_mean_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, v850_mean, lat,lon, time);
    A11_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A11, lat,lon, time);
    A21_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A21, lat,lon, time);
    A22_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A22, lat,lon, time);
    A31_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A31, lat,lon, time);
    A33_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A33, lat,lon, time);
    A42_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A42, lat,lon, time);
    A43_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A43, lat,lon, time);
    A44_t = interp3(wind_xyz_lat, wind_xyz_lon, wind_xyz_time, A44, lat,lon, time);
    F1_t = fourier_series(time);
    F2_t = fourier_series(time);
    F3_t = fourier_series(time);
    F4_t = fourier_series(time);
    %calculate the u v at the time step and location
    u200 = u200_mean_t + A11_t * F1_t;
    u850 = u850_mean_t + A31_t * F1_t + A33_t * F3_t;
                             % 1  2    3             4     5       6    7
    plot_wind (step, :) = [time, u200, u200_mean_t, A11_t, F1_t, A21_t, F2_t, A22_t, A31_t, A33_t, F3_t, A42_t, A43_t, A44_t,F4_t];
    time = time + 0.5/24;
    step = step + 1;
end

%plot u200
plot(plot_wind(:,1), plot_wind(:,2));
hist(plot_wind(:,2));
plot(plot_wind(1:48:end,1), plot_wind(1:48:end,2));
plot(plot_wind(:,1), (plot_wind(:,4).* plot_wind(:,5))./plot_wind(:,3));
hist(plot_wind(:,4).* plot_wind(:,5)./plot_wind(:,3));
hist(plot_wind(:,4)./plot_wind(:,3));
hist(plot_wind(:,4));
plot(plot_wind(1:8,1), plot_wind(1:8,2));
%plot u200_mean
plot(plot_wind(:,1), plot_wind(:,3));
hist(plot_wind(:,3));
%plot A11
plot(plot_wind(:,1), plot_wind(:,4));
%plot F1
plot(plot_wind(:,1), plot_wind(:,5));
hist(plot_wind(:,5));

hist(plot_wind(:,5).*plot_wind(:,4));

%plot u850


%plot example


%% plot all the historical tracks
%plot all the tracks
nc_f = netcdf.open('attracks.nc','NOWRITE');
lat = netcdf.getVar(nc_f, 6, [1,1],[190,1752]);
lon = netcdf.getVar(nc_f, 7, [1,1],[190,1752]);

%covert zero to nan
lat(lat == 0) = NaN;
lon(lon < 200) = NaN;

for i = 1:1752%699
    plot( lon(:,i), lat(:,i),'Color', rand(1,3));
    hold on
end
axis([240 380 5 55])
plot_coast(250:350, 5:50)
grid off

%% wind 

%% plot genesis locations
scatter(sample(:,2),sample(:,1),'filled');
axis([250 350 5 50]);
hold on
plot_coast(250:350, 5:50)

%% plot the standard diviation and mean
pcolor(u200_mean(:,:,8)');
pcolor(A11(:,:,8)');
pcolor(u850_mean(:,:,8)');
pcolor(v850_mean(:,:,8)');
pcolor(A44(:,:,8)');
pcolor(A44(:,:,8)'./v850_mean(:,:,8)');
pcolor(v850_mean(:,:,8)'./u850_mean(:,:,8)');
pcolor(A11(:,:,8)');

caxis([-3 3])

%location plot
scatter(array_record(:,2),array_record(:,1),'filled');
axis([250 350 5 50]);

%% clustering  http://www.datalab.uci.edu/resources/CCT/doc/demo.html   #####################
% read historical data
%plot all the tracks
nc_f = netcdf.open('attracks.nc','NOWRITE');
lat_hist = netcdf.getVar(nc_f, 6, [1,1],[190,1752]);
lon_hist = netcdf.getVar(nc_f, 7, [1,1],[190,1752]);

%covert zero to nan
lat_hist(lat_hist == 0) = NaN;
lon_hist(lon_hist < 200) = NaN;

%define a cell array to store tracks
Y  = cell(size(lat_hist, 2), 1);
for i = 1:size(lat_hist, 2)
    Y{i} = [squeeze(lon_hist(:,i)-360),squeeze(lat_hist(:,i))];
end

% clustering historical data


%read sythetic tracks

% clustering sythetic tracks

%% %% select tracks that pass location of interest in historical tracks 

nc_f = netcdf.open('attracks.nc','NOWRITE');
year = netcdf.getVar(nc_f, 0, 1,1752);
%starting from  index 1035 is tracks after 1970
count  = 1752-1035+1;
lat_hist = netcdf.getVar(nc_f, 6, [1,1035],[190,count]);
lon_hist = netcdf.getVar(nc_f, 7, [1,1035],[190,count]);

%covert zero to nan
lat_hist(lat_hist == 0) = NaN;
lon_hist(lon_hist < 200) = NaN;

tracks_miami_hist = [];
tracks_boston_hist = [];
tracks_norfolk_hist = [];
tracks_newyork_hist = [];
for i = 1: count
    for step = 1:190 
        lat = lat_hist(step,i);
        lon = lon_hist(step,i);
        if lat<42.36+1 & lat>42.36-1 & lon<289.94+1 & lon>289.94-1
            tracks_boston_hist (end+1) = i;
        elseif  lat<25.76+1 & lat>25.76-1 & lon<279.81+1 & lon>279.81-1
            tracks_miami_hist (end+1) = i;
        elseif  lat<36.85+1 & lat>36.85-1 & lon<283.7+1 & lon>283.7-1
            tracks_norfolk_hist (end+1) = i; 
        elseif  lat<40.71+1 & lat>40.71-1 & lon<286+1 & lon>286-1
            tracks_newyork_hist (end+1) = i; 
        end
    end
end

% remove duplicate in the array
tracks_boston_hist = unique (tracks_boston_hist);
tracks_miami_hist = unique (tracks_miami_hist);
tracks_norfolk_hist = unique (tracks_norfolk_hist);
tracks_newyork_hist = unique (tracks_newyork_hist);


%plot coast
plot_coast(250:350, 5:50,[.7 .7 .7])
hold on
grid off;

% plot all the tracks as grey
for i =  1: count
    plot(lon_hist(:,i), lat_hist(:,i), 'Color', [0.8 0.8 0.8] ); %'Color' , rand(1,3)  'Color', [17 17 17]
    hold on
end

%plot the tracks within interest area
% for i =  1: length(tracks_newyork)
%     %plot(sythetic_tracks(tracks_norfolk(i),:,3), sythetic_tracks(tracks_norfolk(i),:,2), 'linewidth', 1.5, 'Color', rand(1,3));
%     %plot(sythetic_tracks(tracks_miami(i),:,3), sythetic_tracks(tracks_miami(i),:,2), 'linewidth', 1.5, 'Color', rand(1,3));
%     %plot(sythetic_tracks(tracks_boston(i),:,3), sythetic_tracks(tracks_boston(i),:,2), 'linewidth', 1.5, 'Color', rand(1,3));
%     plot(sythetic_tracks(tracks_newyork(i),:,3), sythetic_tracks(tracks_newyork(i),:,2), 'linewidth', 1.5, 'Color', rand(1,3));
%     hold on
% end
% set(gcf,'color','w')
% scatter(286,40.71,200,'filled','r')


%norfolk
for i =  1: length(tracks_norfolk_hist)
    plot(lon_hist(:, tracks_norfolk_hist(i)), lat_hist(:,tracks_norfolk_hist(i)), 'linewidth', 1.5, 'Color', rand(1,3));
    hold on
end
set(gcf,'color','w')
scatter(283.7,36.85,200, 'filled','r')

    