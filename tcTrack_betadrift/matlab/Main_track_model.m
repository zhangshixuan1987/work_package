%Propagate
clear variables;
clc;
tic;
%% 1. load gensis
load ('sample_50000');
% how many tracks do we want to run?
N= 50000;
time_step = 0.5/24;  %unit days
%% 2. load synthetic time series of wind

% CMIP6 DATA:
mod_names = {'GFDLCM4', 'CanESM5', 'MIROC6', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'IPSL-CM6A-LR', 'EC-Earth3', 'CMCC-CM2-SR5'};
for imod = 2:8
    model_name = mod_names{imod};
    
    load(['wind_data_', model_name])

    % for boundary cells that are nans, convert them to 0 so that the model
    % does not run into error. The boundaries should not affect model results since they were cutoff:
    beta_drift(isnan(beta_drift)) = 0;

    load('Fourier_series')

    %% 3 create wind time series
    %create array of u_mean, As, Fs. Enlarge them to be 41x21x182. Have a value
    %for each day.
    %Then combine them, to create u_ts and v_ts
    %3d interplolate to get u, v at the time and location

    days_36500 = 1:36500;
    %time_bound = [0 31 59 90 120 151 181 212 243 273 304 334 365];
    %time_bound = [0 15 30 45 60 75 90 106 121 136 152 167 182];
    time_bound = [0 31 59 90 120 151 181 212 243 273 304 334 365] * 36500/365;
    u200_mean_36500 = zeros(41,21,36500);
    v200_mean_36500 = zeros(41,21,36500);
    u850_mean_36500 = zeros(41,21,36500);
    v850_mean_36500 = zeros(41,21,36500);
    A11_36500 = zeros(41,21,36500);
    A21_36500 = zeros(41,21,36500);
    A22_36500 = zeros(41,21,36500);
    A31_36500 = zeros(41,21,36500);
    A33_36500 = zeros(41,21,36500);
    A42_36500 = zeros(41,21,36500);
    A43_36500 = zeros(41,21,36500);
    A44_36500 = zeros(41,21,36500);
    F1_36500 = zeros(41,21,36500);
    F2_36500 = zeros(41,21,36500);
    F3_36500 = zeros(41,21,36500);
    F4_36500 = zeros(41,21,36500);
    wind_xyz_lat_36500 = zeros(41,21,36500);
    wind_xyz_lon_36500 = zeros(41,21,36500);
    wind_xyz_time_36500 = zeros(41,21,36500);

    for i = 1:41
        for j = 1:21
            for k = 1:12
                u200_mean_36500(i,j,(time_bound(k)+1):time_bound(k+1))=u200_mean(i,j,k);
                v200_mean_36500(i,j,(time_bound(k)+1):time_bound(k+1))=v200_mean(i,j,k);
                u850_mean_36500(i,j,(time_bound(k)+1):time_bound(k+1))=u850_mean(i,j,k);
                v850_mean_36500(i,j,(time_bound(k)+1):time_bound(k+1))=v850_mean(i,j,k);
                A11_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A11(i,j,k);
                A21_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A21(i,j,k);
                A22_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A22(i,j,k);
                A31_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A31(i,j,k);
                A33_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A33(i,j,k);
                A42_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A42(i,j,k);
                A43_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A43(i,j,k);
                A44_36500(i,j,(time_bound(k)+1):time_bound(k+1))=A44(i,j,k);
                wind_xyz_lat_36500(i,j,(time_bound(k)+1):time_bound(k+1))=wind_xyz_lat(i,j,k);
                wind_xyz_lon_36500(i,j,(time_bound(k)+1):time_bound(k+1))=wind_xyz_lon(i,j,k);
            end
        end
    end
    %fourier series is only depend on the small t
    for m = 1:36500
        wind_xyz_time_36500(:, :, m)=m/36500*365;
        F1_36500(:,:,m)=F1(mod(m,1500)+1);
        F2_36500(:,:,m)=F2(mod(m,1500)+1);
        F3_36500(:,:,m)=F3(mod(m,1500)+1);
        F4_36500(:,:,m)=F4(mod(m,1500)+1);
    end

    u200_36500 = u200_mean_36500 + A11_36500.*F1_36500;
    v200_36500 = v200_mean_36500 + A21_36500.*F1_36500 + A22_36500.*F2_36500;
    u850_36500 = u850_mean_36500 + A31_36500.*F1_36500 + A33_36500.*F3_36500;
    v850_36500 = v850_mean_36500 + A42_36500.*F2_36500 + A43_36500.*F3_36500 + A44_36500.*F4_36500;

    %% 4. propagation
    %Maximum time is 30 days, equals to 1440 time steps
    max_steps = ceil(30/time_step);
    sythetic_tracks = nan (N,max_steps, 9);
    % Load coastal data
    coast = load('coast.mat');
    [Z, R] = vec2mtx(coast.lat, coast.long, ...
        1, [-90 90], [-180 180], 'filled');
    %follow each track from sample
    parfor i = 1: N    
        %initial condition
        step = 1;
        lat = sample(i,1);
        lon = sample(i,2);
        time = sample(i,3);
        time_on_land = 0;
        temp_sythetic_tracks = nan (1,max_steps, 9);
        %for each time step
        while  lat >4 && lat < 50 && lon >260 && lon <355 && step < max_steps  && time_on_land<4  %(u^2 +v^2)>13^2
            %use 2d interpolation to see whether it saves time
            i_time = ceil(time*36500/365);
            %continue if the track survived to the next year
            if i_time> 36500
                i_time = mod(i_time, 36500);
            end
            %calculate the u v at the time step and location
            u200= interp2(wind_xyz_lat_36500(:,:,i_time), wind_xyz_lon_36500(:,:,i_time),u200_36500(:,:,i_time), lat,lon);
            v200= interp2(wind_xyz_lat_36500(:,:,i_time), wind_xyz_lon_36500(:,:,i_time),v200_36500(:,:,i_time), lat,lon);
            u850= interp2(wind_xyz_lat_36500(:,:,i_time), wind_xyz_lon_36500(:,:,i_time),u850_36500(:,:,i_time), lat,lon);
            v850= interp2(wind_xyz_lat_36500(:,:,i_time), wind_xyz_lon_36500(:,:,i_time),v850_36500(:,:,i_time), lat,lon);
            %calculate beta_u, beta_v
            lat_index = ceil((lat-0)/2.5);
            lon_index = ceil((lon-260)/2.5);
            beta_u = beta_drift(lon_index,lat_index,1);
            beta_v = beta_drift(lon_index,lat_index,2);
            %calculate u, v
            alfa = 0.8;
            u = alfa * u850 + (1-alfa) * u200 + beta_u; %added lat beta
            v = alfa * v850 + (1-alfa) * v200 + beta_v;  %unit m/s %2.5
            %calculate new lat lon and convert from meter to degrees using
            %conversion formula
            lat = lat + v * time_step * 24 * 60 * 60 * degree_per_meter(lat, 'lat'); 
            lon = lon + u * time_step * 24 * 60 * 60 * degree_per_meter(lat, 'lon');
    %         %count how long has the track being on land
    %         if lon>180
    %             lon_test_ocean = lon-360;
    %         else
    %             lon_test_ocean = lon;
    %         end
    %         val = ltln2val(Z, R, lat, lon_test_ocean);
    %         isOcean = val == 2;
    %         if isOcean== 0
    %             time_on_land = time_on_land + time_step;
    %         end
            %save in each track's each step in sythetic_tracks
            temp_sythetic_tracks(1, step, :) = [time, lat, lon, u, v,u200,v200,u850,v850];
            step = step + 1;
            time = time + time_step;
        end 
        sythetic_tracks(i, :, :) = temp_sythetic_tracks(1, :, :);
        fprintf('finished %d out of %d tracks.\n',i,N);  %to plot the number of tracks so that we know where we are
    end
    %plot all the tracks
    plot_coast(260:360, 5:60, [0.1 0.7 0.3]) %dark green backgrond
    grid off
    hold on
    for i =  40:100
        plot(sythetic_tracks(i,:,3), sythetic_tracks(i,:,2), 'Color', rand(1,3) ); %'Color' , rand(1,3)  'Color', [0.8 0.8 0.8]
    end

    print(['random 60 tracks_', model_name], '-dpng')
    close
    toc;

    fig_handle = figure('PaperUnits', 'inches', 'PaperPosition', [0,0,10,4.5]);
    %plot coast
    plot_coast(250:350, 5:50,[.7 .7 .7])
    hold on
    grid off;
    lat_hist = sythetic_tracks(:,:,2);
    lon_hist = sythetic_tracks(:,:,3);
    % plot all the tracks as grey
    for i =  100: 2000
        plot(lon_hist(i,:), lat_hist(i,:), 'Color', [0.8 0.8 0.8] ); %'Color' , rand(1,3)  'Color', [17 17 17]
        hold on
    end
    set(gcf,'color','w')
    print(['plot_first_100-2000_tracks_', model_name],'-dpng')
    close

    %plot the translation vector on top
    load('beta_drift', 'translation_vector')
    for lon = 1:41
        for lat = 1:21
            lat_2d(lon, lat) = 0+(lat-1)*2.5;
            lon_2d(lon, lat) = 260 + (lon-1)*2.5;
        end
    end
    scale_factor = 0.4;
    quiver(lon_2d, lat_2d, squeeze(translation_vector(:,:,1))*scale_factor,...
        squeeze(translation_vector(:,:,2))*scale_factor,'AutoScale', 'off', 'color','b');
    quiver(251, 7, 12*scale_factor,0, 'AutoScale', 'off', 'color','b')
    text(251, 6, '12 m/s')
    print(['plot_tracks_and_historical_tranlation_vector_', model_name],'-dpng')
    close
    
    synthetic_tracks = sythetic_tracks(:,1:12:end, :);
    save (['sythetic_tracks_new_', model_name], 'synthetic_tracks');
    toc;

end

% mod_names = {'GFDLCM4', 'CanESM5', 'MIROC6', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'IPSL-CM6A-LR', 'EC-Earth3', 'CMCC-CM2-SR5'};
% for imod = 1:8
%     model_name = mod_names{imod};
%     disp(model_name);
%     load (['sythetic_tracks_', model_name], 'synthetic_tracks');
%     save (['sythetic_tracks_new_', model_name], 'synthetic_tracks');
% end
