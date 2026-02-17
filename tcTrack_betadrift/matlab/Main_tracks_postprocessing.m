clear variables; close all; clc
tic
whether_plot_each_city = 'n'; %'y', or 'n'

%% 0.0 calculate historical tracks displacement and statistics
% nc_f = netcdf.open('attracks.nc','NOWRITE');
% ncdisp('attracks.nc')
% latmc = ncread('attracks.nc', 'latmc');
% longmc = ncread('attracks.nc', 'longmc');
% yearic = ncread('attracks.nc', 'yearic');
% monthmc = ncread('attracks.nc', 'monthmc');
% daymc = ncread('attracks.nc', 'daymc');
% hourmc = ncread('attracks.nc', 'hourmc');
% 
% %define the array for storage the displacement
% displace_x_obs = nan((1753-1036)*191, 1);
% displace_y_obs = nan((1753-1036)*191, 1);
% count = 1;
% for i = 1036:1753
%     lat_old = latmc(1, i);
%     lon_old = longmc(1, i);
%     for j = 2:191
%         lat = latmc(j, i);
%         lon = longmc(j, i);
%         if lon==0  %If the track ends
%             break
%         elseif lat<10 || lat >30 || lon<280 || lon>330
%             break
%         else
%             displace_x_obs(count) = (lon-lon_old)/degree_per_meter(lat, 'lon')*degree_per_meter(lat, 'lat');
%             displace_y_obs(count) = lat-lat_old;
%             count = count+1;
%             lon_old= lon;
%             lat_old= lat;
%         end
%     end
% end
% 
% %remove of nan values
% idx = find(isnan(displace_x_obs));
% displace_x_obs(idx) = [];
% displace_y_obs(idx) = [];
% 
% %% 0.1 calculate historical tracks' possibility of passing U.S. continent
% disp('0.1 calculate historical tracks possibility of passing U.S. continent')
% %0.1.1 record tracks and their time of passing U,S, continent
% %0.1.1.1 create grid %lat: 5-55, lon -105 (255) - -15(345)
% resolution = 2;
% dim_lat = (55-5)/resolution;
% dim_lon = (345-255)/resolution;
% lat_grid = 5+resolution/2:resolution:55-resolution/2; %the lat and lon of the grid
% lon_grid = 255+resolution/2:resolution:345-resolution/2;
% pass_tracks = cell(dim_lon, dim_lat);% the tracks that pass the grid
% pass_hit_tracks = cell(dim_lon, dim_lat);% the tracks that pass the grid that hit U.S.
% pass_track_count = nan(dim_lon, dim_lat); % count of the tracks that pass the grid
% pass_hit_track_count = nan(dim_lon, dim_lat); % count of the tracks that pass the grid and pass U.S.
% %0.1.1.2 record whether the tracks pass by U.S. and the time
% pass_us = zeros(1753,1);  % Whether a track pass U.S.
% pass_us_time = nan(1753,1); %the time that pass U.S.
% contusa = shaperead('C:\2016\20_hurricane\GIS\ne_10m_admin_0_countries\CONTINENTAL_USA_continous.shp');
% contusa.X = contusa.X + 360;
% 
% for i = 1:1753
%     for j = 1:191
%         lat = latmc(j, i);
%         lon = longmc(j, i);
%         if lat ==0
%             break
%         end
%         if inpolygon(lon, lat, contusa.X, contusa.Y)==1
%             pass_us(i) =1;
%             pass_us_time(i) = datenum(yearic(i), monthmc(j,i), daymc(j,i), hourmc(j,i),0,0);
%             break
%         end
%     end
% end
% 
% %0.1.2 (record the tracks that pass each grid) .
% for i = 1:1753
%     for j = 1:191
%         lat = latmc(j, i);
%         lon = longmc(j, i);
%         %if track endand move to the next track
%         if lon==0 ||lat>55 || lat<5 || lon>345 || lon<255
%             break
%         end
%         latgrid_idx = ceil((lat-5)/resolution);  %calculate the grid index of latitude
%         longrid_idx = ceil((lon-255)/resolution);  %calculate the grid index of longitude
%         pass_tracks{longrid_idx,latgrid_idx}  = [pass_tracks{longrid_idx,latgrid_idx}, i]; %add the tracks to array
%         if pass_us(i) ==1
%             pass_hit_tracks{longrid_idx,latgrid_idx}  = [pass_hit_tracks{longrid_idx,latgrid_idx}, i]; %add the tracks to array
%         end
%         %if track pass us, stop recording and move to the next track
%         if inpolygon(lon, lat, contusa.X, contusa.Y)==1
%             break
%         end
%     end
% end
% 
% %0.1.3 calculate the possbility of passing U.S
% for longrid_idx = 1: dim_lon
%     for latgrid_idx = 1:dim_lat
%         if numel(pass_tracks{longrid_idx,latgrid_idx})>=10
%             pass_track_count(longrid_idx, latgrid_idx) = numel(unique(pass_tracks{longrid_idx,latgrid_idx}));
%             pass_hit_track_count(longrid_idx, latgrid_idx) = numel(unique(pass_hit_tracks{longrid_idx,latgrid_idx}));
%         end
%     end
% end
% 
% [lat_mesh, lon_mesh] = meshgrid(lat_grid, lon_grid);
% pass_hit_risk = pass_hit_track_count./pass_track_count; %the risk of the track will landfall
% 
% %%plot the possibility of a track eventually passing over USA
% fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 10 7]);
% contourf(lon_mesh, lat_mesh, pass_hit_risk*100)
% colormap('jet')
% %colormapeditor
% cbar = colorbar;
% ylabel(cbar, 'Probability(%)')
% cbar.Location = 'southoutside';
% hold on
% plot_coast(255:345, 5:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
% grid off
% plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
% title('Possibility of a track eventually passing over USA')
% print('historical tracks hit possibility', '-dpng')
% close
% 
% %% 0.2 calculate historical tracks' timing of passing U.S. continent
% disp('0.2 calculate historical tracks timing of passing U.S. continent')
% %0.2.1 record tracks and their time of passing U,S, continent
% %0.2.1.1 create grid %lat: 5-55, lon -105 (255) - -15(345)
% resolution_time = 1;
% dim_lat_time = (55-5)/resolution_time;  %the dimension at different resolution_time for timing of passing usa
% dim_lon_time = (345-255)/resolution_time;
% lat_grid_time = 5+resolution_time/2:resolution_time:55-resolution_time/2; %the lat and lon of the grid
% lon_grid_time = 255+resolution_time/2:resolution_time:345-resolution_time/2;
% pass_hit_timelength = cell(dim_lon_time, dim_lat_time);% the tracks that pass the grid that hit U.S.
% pass_hit_timelength_mean = nan(dim_lon_time, dim_lat_time); %the average lenth of time that track land from this position
% 
% %0.2.2 (record the tracks that pass each grid) .
% for i = 1:1753
%     for j = 1:191
%         lat = latmc(j, i);
%         lon = longmc(j, i);
%         %if track endand move to the next track
%         if lon==0 ||lat>55 || lat<5 || lon>345 || lon<255
%             break
%         end
%         latgrid_idx = ceil((lat-5)/resolution_time);  %calculate the grid index of latitude
%         longrid_idx = ceil((lon-255)/resolution_time);  %calculate the grid index of longitude
%         if pass_us(i) ==1
%             timenow = datenum(yearic(i), monthmc(j,i), daymc(j,i), hourmc(j,i),0,0);
%             pass_hit_timelength{longrid_idx,latgrid_idx} = [pass_hit_timelength{longrid_idx,latgrid_idx}, pass_us_time(i)-timenow];
%         end
%         %if track pass us, stop recording and move to the next track
%         if inpolygon(lon, lat, contusa.X, contusa.Y)==1
%             break
%         end
%     end
% end
% 
% %0.2.3 calculate the possbility of passing U.S, and the average timing of passing U.S.
% for longrid_idx = 1: dim_lon_time
%     for latgrid_idx = 1:dim_lat_time
%         if numel(pass_hit_timelength{longrid_idx,latgrid_idx})>=5
%              pass_hit_timelength_mean(longrid_idx, latgrid_idx) = nanmean(pass_hit_timelength{longrid_idx,latgrid_idx});
%         end
%     end
% end
% 
% [lat_mesh_time, lon_mesh_time] = meshgrid(lat_grid_time, lon_grid_time);
% pass_hit_risk = pass_hit_track_count./pass_track_count*100; %percentage the risk of the track will landfall
% 
% %%Plot the average length of time a tropical cyclone takes to make landfall
% %%at USA
% pass_hit_timelength_mean_below14 = pass_hit_timelength_mean;
% pass_hit_timelength_mean_below14(find(pass_hit_timelength_mean>14)) =nan;
% 
% fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 10 7]);
% %contourf(lon_mesh_time, lat_mesh_time, pass_hit_timelength_mean_below14)
% pplot = pcolor(lon_mesh_time, lat_mesh_time, pass_hit_timelength_mean_below14);
% set(pplot, 'EdgeColor', 'none');
% colormap('hsv')
% %colormapeditor
% %colormap(flipud(colormap))
% cbar = colorbar;
% ylabel(cbar, 'days')
% cbar.Location = 'southoutside';
% hold on
% plot_coast(255:345, 5:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
% grid off
% plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
% title('Average length of time a tropical cyclone takes to make landfall')
% % xlabel('longitude')
% % ylabel('latitude')
% print('historical tracks length of time before hit', '-dpng')
% close
% 
% save('historical_tracks_statistics_saved')

load('historical_tracks_statistics_saved')
%% 0.4 load historical tracks that pass by the 52 cities
%save('historical_hit_52_cities', 'city_names', 'city_lat', 'city_lon', 'city_hit_count_obs', 'city_hit_risk_obs', 'city_tracks_obs');
load('historical_hit_52_cities')

%% 1.0 load synthetic tracks
disp('1.0 load synthetic tracks')
load('sythetic_tracks')  %[time, lat, lon, u, v, time_on_land];
N = size(sythetic_tracks, 1);

%% 1.1 compare the 52 city hit risks
disp('1.1 compare the 52 city hit risks')
city_hit_risk_obs = city_hit_risk_obs*100;
city_hit_risk = city_hit_count / N *100;
%fit linear regression
p = polyfit(city_hit_risk, city_hit_risk_obs,1);
risk_fit = polyval(p, city_hit_risk);
risk_resid = city_hit_risk_obs - risk_fit;
SSresid = sum(risk_resid.^2);
SStotal = (length(city_hit_risk_obs)-1)*var(city_hit_risk_obs);
rsq = 1 - SSresid/SStotal;

fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 4 4]);
scatter(city_hit_risk, city_hit_risk_obs)
hold on
plot(city_hit_risk, risk_fit, 'LineWidth', 1)
text(3.2, 1.5, ['y = x * ' num2str(p(1)) '+' num2str(p(2))])
text(3.2, 1.2, ['R square = ' num2str(rsq)])
ylabel('Observed risk(%)')
xlabel('Modeled risk(%)')
ylim([0 5])
xlim([0 5])
print('city_hit_risk', '-dpng')
close

if strcmp(whether_plot_each_city, 'y')
    disp('    plot the 52 city hit tracks, modeled vs observed...')
    for i = 1: length(city_hit_count)
        figname = ['More_Plots\city_' num2str(i) '_' city_names{i} '_mod vs obs tracks' ];
        fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 7 3.5]);
        plot_coast(255:345, 5:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
        grid off
        hold on
        %plot modeled tracks in grey
        for j = 1:numel(city_tracks{i})
            plot(sythetic_tracks(city_tracks{i}(j), :, 3),sythetic_tracks(city_tracks{i}(j), :, 2), 'Color', [0.8 0.8 0.8])
        end
        %plot observed tracks in blue
        latmc_plot = latmc;
        longmc_plot = longmc;
        latmc_plot(find(latmc ==0)) = nan;
        longmc_plot(find(longmc ==0)) = nan;
        for j = 1: numel(city_tracks_obs{i})
            plot(longmc_plot(:, city_tracks_obs{i}(j)),latmc_plot(:, city_tracks_obs{i}(j)), 'b')
        end
        %plot the city as red dot
        plot(city_lon(i), city_lat(i), 'ro')
        text(320, 10, city_names{i})
        print(figname, '-dpng')
        close
    end
end

%% 1.2 calculate displacement at 6 hr interval & RECORD whether track pass US
disp('1.2 load synthetic tracks & calculate displacement at 6 hr interval & RECORD whether track pass US')
time_step = 0.5/24;  %unit days
sixhr_steps = 6/24/time_step;
max_steps = size(sythetic_tracks, 2);
%define the array for storage the displacement
displace_x = nan(N*max_steps, 1);
displace_y = nan(N*max_steps, 1);

% define the array to save possibilities and track count 
pass_tracks_mod = cell(dim_lon, dim_lat);% the tracks that pass the grid
pass_hit_tracks_mod = cell(dim_lon, dim_lat);% the tracks that pass the grid that hit U.S.
pass_track_count_mod = nan(dim_lon, dim_lat); % count of the tracks that pass the grid
pass_hit_track_count_mod = nan(dim_lon, dim_lat); % count of the tracks that pass the grid and pass U.S.
%0.1.1.2 record whether the tracks pass by U.S. and the time
pass_us_mod = zeros(N,1);  % Whether a track pass U.S.
pass_us_time_mod = nan(N,1); %the time that pass U.S.
pass_us_step_mod = nan(N,1); %the time step that pass U.S.
contusa = shaperead('C:\2016\20_hurricane\GIS\ne_10m_admin_0_countries\CONTINENTAL_USA_continous.shp');
contusa.X = contusa.X + 360;

%define the array to save the passing time length at 1 degree resolution
pass_hit_timelength_mod = cell(dim_lon_time, dim_lat_time);% the tracks that pass the grid that hit U.S.
pass_hit_timelength_mean_mod = nan(dim_lon_time, dim_lat_time); %the average lenth of time that track land from this position


count = 1;
for i = 1:N
    %record whether synthetic tracks pass U.S.
    lats = sythetic_tracks(i, 6:6:1440, 2);%6 steps each time to save time
    lons = sythetic_tracks(i, 6:6:1440, 3);
    idx_inusa = inpolygon(lons, lats, contusa.X, contusa.Y); 
    if nansum(idx_inusa)>0  %if there are us landfalls
        pass_us_mod(i) =1;
        idx_inusa_list = find(idx_inusa==1);
        pass_us_time_mod(i) = sythetic_tracks(i, idx_inusa_list(1)*6, 1);  %6 steps each time to save time
        pass_us_step_mod(i) = idx_inusa_list(1)*6;
    end
    
    %calculate displacement
    lat_old = sythetic_tracks(i, 1, 2);
    lon_old = sythetic_tracks(i, 1, 3);
    for j = (sixhr_steps+1):sixhr_steps:max_steps
        lat = sythetic_tracks(i, j, 2);
        lon = sythetic_tracks(i, j, 3);
        if isnan(lat)
            break
        elseif lat<10 || lat >30 || lon<280 || lon>330
            break
        else
            displace_x(count) = (lon-lon_old)/degree_per_meter(lat, 'lon')*degree_per_meter(lat, 'lat');
            displace_y(count) = lat-lat_old;
            count = count+1;
            lon_old= lon;
            lat_old= lat;
        end
    end
end

%remove of nan values
idx = find(isnan(displace_x));
displace_x(idx) = [];
displace_y(idx) = [];
toc
%% 1.3 PLOT compare displacement histgram
disp('1.3 PLOT compare displacement histgram')
fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 6 6]);
%%plot zonal displacement
% randomly select the same number of tracks as observation to plot
% histogram
subplot(2,1,1)
num_comparison = size(displace_x_obs, 1);
num_mod_displ = size(displace_x, 1);
idx= randi(num_mod_displ, [num_comparison,1]);
nbins=40;
[hist_obs,centers] = hist(displace_x_obs,nbins);
[hist_mod] = hist(displace_x(idx),centers);
width1 = 0.35;
bar(centers,hist_mod,width1,'FaceColor','r',....
                     'EdgeColor','none');
hold on
width2 = width1;
bar(centers+0.05,hist_obs,width2,'FaceColor','b',...
                     'EdgeColor','b');
text(-2.8,550, '(a)')
hold off
xlim([-3,3])
xlabel('6 hour zonal displacements(deg lat)')
ylabel('Number of occurances')
legend('Model','Best track since 1970') % add legend
disp(['the mean of modeled displacement in x direction is ' num2str(mean(displace_x(idx))) ', while observed is ' num2str(mean(displace_x_obs))])

%%plot meridional displacement
% randomly select the same number of tracks as observation to plot
% histogram
subplot(2,1,2)
[hist_obs,centers] = hist(displace_y_obs,nbins);
[hist_mod] = hist(displace_y(idx),centers);
bar(centers,hist_mod,width1,'FaceColor','r',....
                     'EdgeColor','none');
hold on
width2 = width1;
bar(centers+0.05,hist_obs,width2,'FaceColor','b',...
                     'EdgeColor','b');
text(-1.8,900, '(b)')
hold off
xlim([-2,3])
xlabel('6 hour meridional displacements(deg lat)')
ylabel('Number of occurances')
legend('Model','Best track since 1970') % add legend
print('histogram_zonal and meridional_displacements','-dpng')
close
disp(['the mean of modeled displacement in y direction is ' num2str(mean(displace_y(idx))) ', while observed is ' num2str(mean(displace_y_obs))])
toc
%% 1.4 record the probability of TC passing grid, and the time of passing
disp('1.4 record the probability of TC passing grid, and the time of passing ')
tic
%POSSIBILITY OF modeled track passing by 
for i = 1:N
    for j = 1:max_steps
        lat = sythetic_tracks(i, j, 2);
        lon = sythetic_tracks(i, j, 3);
        if isnan(lat) ||lat>55 || lat<5 || lon>345 || lon<255
            break
        end
        %use the 2 degree resolution to record the passing tracks
        latgrid_idx = ceil((lat-5)/resolution);  %calculate the grid index of latitude
        longrid_idx = ceil((lon-255)/resolution);  %calculate the grid index of longitude
        pass_tracks_mod{longrid_idx,latgrid_idx}  = [pass_tracks_mod{longrid_idx,latgrid_idx}, i]; %add the tracks to array
        if pass_us_mod(i) ==1
            pass_hit_tracks_mod{longrid_idx,latgrid_idx}  = [pass_hit_tracks_mod{longrid_idx,latgrid_idx}, i]; %add the tracks to array
        end
        %use the 1 degree resolution to record the passing time length
        latgrid_idx = ceil((lat-5)/resolution_time);  %calculate the grid index of latitude
        longrid_idx = ceil((lon-255)/resolution_time);  %calculate the grid index of longitude
        if pass_us_mod(i) ==1
            timenow = sythetic_tracks(i, j, 1);
            pass_hit_timelength_mod{longrid_idx,latgrid_idx} = [pass_hit_timelength_mod{longrid_idx,latgrid_idx}, pass_us_time_mod(i)-timenow];
            %when the track passes U.S, move to calculate the next track
            if j == pass_us_step_mod(i)
                break
            end
        end 
    end
end


%calculate the possbility of passing U.S
for longrid_idx = 1: dim_lon
    for latgrid_idx = 1:dim_lat
        if numel(pass_tracks_mod{longrid_idx,latgrid_idx})>=10
            pass_track_count_mod(longrid_idx, latgrid_idx) = numel(unique(pass_tracks_mod{longrid_idx,latgrid_idx}));
            pass_hit_track_count_mod(longrid_idx, latgrid_idx) = numel(unique(pass_hit_tracks_mod{longrid_idx,latgrid_idx}));
        end
    end
end
pass_hit_risk_mod = pass_hit_track_count_mod./pass_track_count_mod*100; %percentage the risk of the track will landfall


%calculate the average timing of passing U.S.
for longrid_idx = 1: dim_lon_time
    for latgrid_idx = 1:dim_lat_time
        if numel(pass_hit_timelength_mod{longrid_idx,latgrid_idx})>=5
             pass_hit_timelength_mean_mod(longrid_idx, latgrid_idx) = nanmean(pass_hit_timelength_mod{longrid_idx,latgrid_idx});
        end
    end
end
toc
%% 1.5 plot the possibility of making landfall at USA
fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 14 4.9]);
subplot(1,2,1)
contourf(lon_mesh, lat_mesh, pass_hit_risk_mod)
colormap('jet')
%colormapeditor
cbar = colorbar;
ylabel(cbar, 'Probability(%)')
cbar.Location = 'southoutside';
hold on
plot_coast(255:345, 5:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
grid off
plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
title('Modeled possibility of a track eventually passing over USA')

subplot(1,2,2)
contourf(lon_mesh, lat_mesh, pass_hit_risk)
colormap('jet')
%colormapeditor
cbar = colorbar;
ylabel(cbar, 'Probability(%)')
cbar.Location = 'southoutside';
hold on
plot_coast(255:345, 5:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
grid off
plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
title('Observed possibility of a track eventually passing over USA')
print('historical tracks hit possibility_compare with model', '-dpng')
close
%% 1.6 plot average length of time TC to make landfall at USA
pass_hit_timelength_mean_below14_mod = pass_hit_timelength_mean_mod;
pass_hit_timelength_mean_below14_mod(find(pass_hit_timelength_mean_mod>14)) =nan;

fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 14 4.9]);
subplot(1,2,1)
pplot = pcolor(lon_mesh_time, lat_mesh_time, pass_hit_timelength_mean_below14_mod);
set(pplot, 'EdgeColor', 'none');
colormap('hsv')
%colormapeditor
%colormap(flipud(colormap))
cbar = colorbar;
ylabel(cbar, 'days')
cbar.Location = 'southoutside';
hold on
plot_coast(255:345, 5:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
grid off
plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
title('Modeled average length of time a tropical cyclone takes to make landfall')

subplot(1,2,2)
pplot = pcolor(lon_mesh_time, lat_mesh_time, pass_hit_timelength_mean_below14);
set(pplot, 'EdgeColor', 'none');
colormap('hsv')
%colormapeditor
%colormap(flipud(colormap))
cbar = colorbar;
ylabel(cbar, 'days')
cbar.Location = 'southoutside';
hold on
plot_coast(255:345, 5:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
grid off
plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
title('Observed length of time a tropical cyclone takes to make landfall')
% xlabel('longitude')
% ylabel('latitude')
print('historical tracks length of time before hit compare with model', '-dpng')
close

%% 1.7 plot the possibility of hurricane hit risk for cities

fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0 0 10 4]);
subplot(1,2,1)
plot_coast(255:300, 20:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
grid off
hold on
plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
scatter(city_lon, city_lat, 20, city_hit_risk, 'filled')
colormap('jet')
caxis([1 4.5]);
cbar = colorbar;
ylabel(cbar, 'risk (%)')
cbar.Location = 'southoutside';
title('Modeled tropical cyclone hit risk')

subplot(1,2,2)
plot_coast(255:300, 20:55, [0.7 0.7 0.7]) %[0.8 0.5 0.1]
grid off
hold on
plot(contusa.X, contusa.Y, 'k', 'LineWidth', 1)
scatter(city_lon, city_lat, 20, city_hit_risk_obs, 'filled')
colormap('jet')
caxis([1 4.5]);
cbar = colorbar;
ylabel(cbar, 'risk (%)')
cbar.Location = 'southoutside';
title('Observed tropical cyclone hit risk')
print('historical city hit risk compare with model', '-dpng')
close

%% save the data
save('post_processing_data', '-v7.3')
toc