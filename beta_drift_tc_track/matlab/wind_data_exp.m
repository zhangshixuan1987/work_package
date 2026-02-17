clear variables;
clc;

% for CMIP6 data, the processed monthly data format is very odd. Models 
% are indexed differently in hist and future files, and some models are
% provided seperately in individual files.
mod_names = {'GFDLCM4', 'CanESM5', 'MIROC6', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'IPSL-CM6A-LR', 'EC-Earth3', 'CMCC-CM2-SR5'};
mod_idx_hist = {6,  13, 10, 11, 12,  9,  1, 15};
mod_idx_futu = {4,   1, 15,  5,   8, 11, 14, 13};
for imod = 1:8
    model_name = mod_names{imod};
    idx_in_file = string(mod_idx_hist{imod});
    UA_name = char('UA'+idx_in_file);
    VA_name = char('VA'+idx_in_file);
    %% read the ASO wind matrix
    if imod<=6
        nfile = 'C:\2016\85_RAFT\CMIP6\ua_200_cmip6_hist.nc';
        u200_aso_raw = squeeze(ncread(nfile,UA_name));
        lats_aso = ncread(nfile, 'LAT73_112');
        lons_aso = ncread(nfile, 'LON139_192');
        nfile = 'C:\2016\85_RAFT\CMIP6\ua_850_cmip6_hist.nc';
        u850_aso_raw = squeeze(ncread(nfile, UA_name));
        nfile = 'C:\2016\85_RAFT\CMIP6\va_200_cmip6_hist.nc';
        v200_aso_raw = squeeze(ncread(nfile, VA_name));
        nfile = 'C:\2016\85_RAFT\CMIP6\va_850_cmip6_hist.nc';
        v850_aso_raw = squeeze(ncread(nfile, VA_name));
    elseif imod == 7
        nfile = 'C:\2016\85_RAFT\CMIP6\ua200_Amon_EC-Earth3_historical_r1i1p1f1_gr_198001-201412_reg.nc';
        u200_aso_raw = squeeze(ncread(nfile,'UA_REG'));
        lats_aso = ncread(nfile, 'LAT73_112');
        lons_aso = ncread(nfile, 'LON139_192');
        nfile = 'C:\2016\85_RAFT\CMIP6\ua850_Amon_EC-Earth3_historical_r1i1p1f1_gr_198001-201412_reg.nc';
        u850_aso_raw = squeeze(ncread(nfile, 'UA_REG'));
        nfile = 'C:\2016\85_RAFT\CMIP6\va200_Amon_EC-Earth3_historical_r1i1p1f1_gr_198001-201412_reg.nc';
        v200_aso_raw = squeeze(ncread(nfile, 'VA_REG'));
        nfile = 'C:\2016\85_RAFT\CMIP6\va850_Amon_EC-Earth3_historical_r1i1p1f1_gr_198001-201412_reg.nc';
        v850_aso_raw = squeeze(ncread(nfile, 'VA_REG'));
    elseif imod == 8
        nfile = 'C:\2016\85_RAFT\CMIP6\ua_Amon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_198001-201412_reg.nc';
        levels = ncread(nfile, 'PLEV');
        ua_aso_raw = ncread(nfile,'UA_REG');
        u200_aso_raw = squeeze(ua_aso_raw(:,:,10,:));
        u850_aso_raw = squeeze(ua_aso_raw(:,:,17,:));
        lats_aso = ncread(nfile, 'LAT73_112');
        lons_aso = ncread(nfile, 'LON139_192');
        nfile = 'C:\2016\85_RAFT\CMIP6\va_Amon_CMCC-CM2-SR5_historical_r1i1p1f1_gn_198001-201412_reg.nc';
        va_aso_raw = ncread(nfile,'VA_REG');
        v200_aso_raw = squeeze(va_aso_raw(:,:,10,:));
        v850_aso_raw = squeeze(va_aso_raw(:,:,17,:));
    end
    %filter extreme values, convert to nan
    u200_aso_raw(u200_aso_raw>100 | u200_aso_raw<-100)=nan;
    v200_aso_raw(v200_aso_raw>100 | v200_aso_raw<-100)=nan;
    u850_aso_raw(u850_aso_raw>100 | u850_aso_raw<-100)=nan;
    v850_aso_raw(v850_aso_raw>100 | v850_aso_raw<-100)=nan;

    % construct 2-d arrays
    [lats_aso_2d, lons_aso_2d] = meshgrid(lats_aso, lons_aso);
    load('..\..\beta_drift\ncept_beta_regression\wind_data', 'lat1', 'lon1')
    [lat_2d, lon_2d] = meshgrid(lat1+2.5/2, lon1+2.5/2);  %2-d array for the ultimate wind input data format

    %interpolate and extrapolate to remove nan values
    for k = 1:size(v850_aso_raw,3)
        xnan = double(u200_aso_raw(:,:,k));
        inan = isnan(xnan);
        zk = griddata(lons_aso_2d(~inan),lats_aso_2d(~inan),xnan(~inan),lons_aso_2d(inan),lats_aso_2d(inan));
        xnan(inan) = zk;
        disp(zk)
        u200_aso_raw(:,:,k) = xnan;
    end

    for k = 1:size(v850_aso_raw,3)
        xnan = double(v200_aso_raw(:,:,k));
        inan = isnan(xnan);
        zk = griddata(lons_aso_2d(~inan),lats_aso_2d(~inan),xnan(~inan),lons_aso_2d(inan),lats_aso_2d(inan));
        xnan(inan) = zk;
        disp(zk)
        v200_aso_raw(:,:,k) = xnan;
    end
    for k = 1:size(v850_aso_raw,3)
        xnan = double(u850_aso_raw(:,:,k));
        inan = isnan(xnan);
        zk = griddata(lons_aso_2d(~inan),lats_aso_2d(~inan),xnan(~inan),lons_aso_2d(inan),lats_aso_2d(inan));
        xnan(inan) = zk;
        disp(zk)
        u850_aso_raw(:,:,k) = xnan;
    end
    for k = 1:size(v850_aso_raw,3)
        xnan = double(v850_aso_raw(:,:,k));
        inan = isnan(xnan);
        zk = griddata(lons_aso_2d(~inan),lats_aso_2d(~inan),xnan(~inan),lons_aso_2d(inan),lats_aso_2d(inan));
        xnan(inan) = zk;
        disp(zk)
        v850_aso_raw(:,:,k) = xnan;
    end
    %% interpolate

    % interpolate to the target format of wind: 
    for k = 1:size(u850_aso_raw,3)
        u850_aso(:,:,k) = interp2( lats_aso_2d,lons_aso_2d,u850_aso_raw(:,:,k), lat_2d, lon_2d);
        u200_aso(:,:,k) = interp2( lats_aso_2d,lons_aso_2d,u200_aso_raw(:,:,k), lat_2d, lon_2d);
        v850_aso(:,:,k) = interp2( lats_aso_2d,lons_aso_2d,v850_aso_raw(:,:,k), lat_2d, lon_2d);
        v200_aso(:,:,k) = interp2( lats_aso_2d,lons_aso_2d,v200_aso_raw(:,:,k), lat_2d, lon_2d);
    end
    % for the longitude at 360, fill the nan values using the last row (this
    % data is not used in the track model because cut off at 355.)
    % ARTIFACTS about CMIP6 DATA, need to repeat the edge since
    % interpolation did not extrapolate
    if isnan(u850_aso(end,2,2))
        u850_aso(end,:,:) = u850_aso(end-1,:,:);
        u200_aso(end,:,:) = u200_aso(end-1,:,:);
        v850_aso(end,:,:) = v850_aso(end-1,:,:);
        v200_aso(end,:,:) = v200_aso(end-1,:,:);
    end
    if isnan(u850_aso(2,end,2))
        u850_aso(:,end,:) = u850_aso(:,end-1,:);
        u200_aso(:,end,:) = u200_aso(:,end-1,:);
        v850_aso(:,end,:) = v850_aso(:,end-1,:);
        v200_aso(:,end,:) = v200_aso(:,end-1,:);
    end
    clear ua va u850_aso_raw u200_aso_raw v850_aso_raw v200_aso_raw
    %% save U200 V200 U850 V850 matrix
    months = 1:12;
    months = months';

    %This is a function that returns a 42by 21 by 12 array of wind speed
    %(mean, and variance, within each of 12 month)
    [u200_mean,u200_var, u200_40_year] = get_wind_matrix( u200_aso);

    %do the same for v200, u850, v850
    [v200_mean, v200_var, v200_40_year] = get_wind_matrix( v200_aso );

    [u850_mean, u850_var, u850_40_year] = get_wind_matrix( u850_aso );

    [v850_mean, v850_var, v850_40_year] = get_wind_matrix( v850_aso );

    %% calculate the covariance
    %returns  3 dimensional matrix
    % v200_40-year is a data set from 1971 jan to 2015 dec 46 years of data
    [A11, A21] = wind_covariance (u200_40_year, v200_40_year);
    [A33, A31] = wind_covariance (u850_40_year, u200_40_year);
    [A22, A42] = wind_covariance (v200_40_year, v850_40_year);
    [A44, A43] = wind_covariance (v850_40_year, u850_40_year);

    %% record the lat, long, time
    wind_xyz_lat = zeros(41,21,12);
    wind_xyz_lon = zeros(41,21,12);
    wind_xyz_time = zeros(41,21,12);

    end_of_month = eomday(2001,1:12)';
    middle_of_month = [cumsum(end_of_month) - end_of_month/2];

    for ts = 1:12
        for lon = 1:41
            for lat = 1:21
                wind_xyz_lat(lon, lat, ts) = 0 + (lat-1)*2.5;
                wind_xyz_lon(lon, lat, ts) = 260 + (lon-1)*2.5;
                wind_xyz_time(lon, lat, ts) = middle_of_month(ts);
            end
        end
    end

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

    beta_drift(:,:,1) =2.567 * dv_dx +...
          0.9285 * f +...
          0.0873 * du_dp +...
          0.2107 * dv_dp +...
         -0.365  * rela_vort_dy +...
         -0.0195 * lon_for_reg +...
         -0.5766 * u850_for_reg +...
         -0.6692 * v850_for_reg +...
         -0.4534;

    beta_drift(:,:,2) =-0.0998 * du_dp +...
          0.2493 * dv_dp +...
          0.1749 * u850_for_reg +...
         -1.2009 * v850_for_reg +...
          1.2707;

    %% save time series wind matrix
    save(['wind_data_', model_name])
end