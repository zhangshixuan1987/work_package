clear variables;
clc;
nc_f = netcdf.open('attracks.nc','NOWRITE');

%[1]display header of the netcdf
ncdisp('attracks.nc')
[ndim,nvar,natt,unlim]=netcdf.inq(nc_f)
[dimname, dimlength] = netcdf.inqDim(nc_f,1)% 1 time length 191
[dimname, dimlength] = netcdf.inqDim(nc_f,0)% 1 ns length 1753
time = netcdf.getVar(nc_f, 1, 0, 191);
[dimname, dimlength] = netcdf.inqDim(nc_f,0)% 1 
[varname, xtype, dimid, natt] = netcdf.inqVar(nc_f,0) %year
[varname, xtype, dimid, natt] = netcdf.inqVar(nc_f,3) %month
[varname, xtype, dimid, natt] = netcdf.inqVar(nc_f,4) %day
[varname, xtype, dimid, natt] = netcdf.inqVar(nc_f,5) %hour
[varname, xtype, dimid, natt] = netcdf.inqVar(nc_f,6) %latmc
[varname, xtype, dimid, natt] = netcdf.inqVar(nc_f,7) %longmc
[varname, xtype, dimid, natt] = netcdf.inqVar(nc_f,8) %vsmc
%netcdf.getVar

%loop that find the first vs >12
%read in all the track records, and save the track genesis locations in Array
array_record = zeros(699,8);
count = 1;
for ns = 0:1752
    for time = 0:190
        velo = netcdf.getVar(nc_f, 8, [time,ns],[1,1]);
        if velo >=12
            year = netcdf.getVar(nc_f, 0, ns,1);
            if year > 1970
                lat = netcdf.getVar(nc_f, 6, [time,ns],[1,1]);
                lon = netcdf.getVar(nc_f, 7, [time,ns],[1,1]);
                month = netcdf.getVar(nc_f, 3, [time,ns],[1,1]);
                day = netcdf.getVar(nc_f, 4, [time,ns],[1,1]);
                hour = netcdf.getVar(nc_f, 5, [time,ns],[1,1]);
                [doy,fraction] = date2doy(datenum(strcat(num2str(month),'/',num2str(day),'/',num2str(year))));
                array_record(count,:) = [lat,lon,year,month,day,hour,doy,velo];   
                count = count + 1;
            end  
            break;    

        end
    end
end

%location plot
% scatter(array_record(:,2),array_record(:,1),'filled');
% axis([250 350 5 50]);

%save a 3d array of lat long, time for 3-d guassian kernel smoothing
  size_lon = (350 - 250)/0.5;  %200
  size_lat = 60/0.5;           %120
  size_days = 365 / 5+1;       %74
  array_3d_genesis_count = zeros(size_lat, size_lon,size_days);
  array_3d_genesis_lat = zeros(size_lat, size_lon,size_days);
  array_3d_genesis_lon = zeros(size_lat, size_lon,size_days);
  array_3d_genesis_day = zeros(size_lat, size_lon,size_days);
for ns = 1:699
    lat = array_record(ns,1);
    lon = array_record(ns,2);
    lat_index = ceil((lat-0)/0.5);
    lon_index = ceil((lon-250)/0.5);
    time_index = ceil(array_record(ns,7)/5);
    array_3d_genesis_count(lat_index,lon_index, time_index) = array_3d_genesis_count(lat_index,lon_index, time_index)+1;

end 
%record lat, long, time values in the big array for future use
for lon_index = 1:size_lon
    for lat_index = 1:size_lat
        for time_index = 1:size_days
            array_3d_genesis_lat(lat_index,lon_index, time_index) = lat_index * 0.5;
            array_3d_genesis_lon(lat_index,lon_index, time_index) = lon_index * 0.5 +250;
            array_3d_genesis_day(lat_index,lon_index, time_index) = time_index * 5;
        end
    end
end
%imagesc([array_3d_genesis_count(:,:,37)])

%% plot the distributition of genesis point every time period and save the figures
%first sort the 
count_historical = zeros(74,1);
for time = 1:74
    plot_coast(250:350, 5:50)
    set(gca,'Ydir','Normal')
    grid off
    hold on
    % go over each track, and plot an dot if it is inside the time period
    for ns = 1:699
        if array_record(ns,7) <time*5 & array_record(ns,7) > time*5-5
            count_historical(time) = count_historical(time) + 1;
            historical_period_track{time}(count_historical(time),1)= array_record(ns,1);
            historical_period_track{time}(count_historical(time),2)= array_record(ns,2);
            scatter(array_record(ns,2), array_record(ns,1));
        end
    end
    title(['time period ', num2str(time), ' number of tracks: ', num2str(count_historical(time))]);
    goodplot(30,20,0.1,11);
    hold off
    pngname = ['More_Plots/time period ', num2str(time), '.png'];
    print(gcf, '-dpng', sprintf('-r%d',120),pngname);
    close;
end

%% calculate possiblity distribution using variable bandwidth guassian kernel
%start the loop go over each location, and count the probablility in each
%box
%set gauss_sum to be a bigger array so that it can take the boundary boxes'
%values. kernel is the final pdf
%maximum of kernel size: 15 degree in lat and lon, which is 15/0.5 = 30
gauss_sum = zeros(200,280,100);
kernel = zeros(120,200,74);

% %calculate sd of 699's location
% sigma_lat_all = std(array_record(:,1))
% sigma_lon_all = std(array_record(:,2))

%follow each track and record the gauss kernel 
for ns = 1: size(array_record,1)
    %load parameters
    lat = array_record(ns,1);
    lon = array_record(ns,2);
    lat_index = ceil((lat-0)/0.5);
    lon_index = ceil((lon-250)/0.5);
    time_index = ceil(array_record(ns,7)/5);
    %default values
    sigma = 5;
    %kernel size
    ksize = 10;
    ksize_time = 1;
    
    %cases when there are less than n = 5 genesis
    if count_historical(time_index) > 1 && count_historical(time_index) < 5
        genesis_subset = historical_period_track{time_index};
        sigma_lat = std(genesis_subset(:,1));
        sigma_lon = std(genesis_subset(:,2));
        sigma = 0.5*(sigma_lat+sigma_lon);
        %silverman approximation
        ksize = 1.06*sigma*size(genesis_subset,1)^(-0.2)*2;
        if ksize>30
            ksize = 30;
        end
    elseif count_historical(time_index) >= 5 
        %find nearest five points
        genesis_subset = historical_period_track{time_index};
        %use knnsearch to find the 5 nearest points to the interest point
        index = knnsearch(genesis_subset,[lat lon],'K', 5);
        %update the subset genesis to only 5 points
        genesis_subset = genesis_subset(index, :);
        sigma_lat = std(genesis_subset(:,1));
        sigma_lon = std(genesis_subset(:,2));
        sigma = 0.5*(sigma_lat+sigma_lon);
        %silverman approximation
        ksize = 1.06*sigma*size(genesis_subset,1)^(-0.2)*2;
        if ksize>30
            ksize = 30;
        end
    end
    %round the ksize to nearest interge
    ksize = ceil(ksize);
    %create a 3 dimensional grid to calculate the distance, and count
    %number of events
    [xg, yg, zg] = ndgrid ( -ksize : ksize, -ksize : ksize, -ksize_time : ksize_time );

    %making the gaussian kernel
    dist = sqrt ( xg .^ 2 + yg .^ 2 + zg .^ 2 ) ;
    gauss = 1 / ( sigma * sqrt ( 2 * pi ) ) * exp ( -0.5 * dist .^ 2 / ( sigma ) ^2 ) ;  
    gauss = gauss / sum ( gauss ( : ) ) ;
    
    gauss_sum(lat_index+30-ksize:lat_index+ksize+30,lon_index+30-ksize:lon_index+ksize+30,time_index:time_index+2*ksize_time)...
        =gauss_sum(lat_index+30-ksize:lat_index+ksize+30,lon_index+30-ksize:lon_index+ksize+30,time_index:time_index+2*ksize_time)...
        +gauss;
    disp(ns);
    disp(ksize);
    disp(sigma);
    disp(sum(gauss(:)));
end

%% remove genesis possibility on land;

% kernel equals to the central part of gauss_sum
kernel = gauss_sum(31:150, 31:230, 2:75);
disp(sum(kernel(:)));

%Determine whether it is land or water
for lat_index = 1:120
    for lon_index = 1:200
        lat = 0.0 + (lat_index-0.5)*0.5;
        lon = 250.0 + (lon_index-0.5)*0.5;
        if landmask (lat, lon-360) == 1
            kernel (lat_index, lon_index,:) = 0;
        end 
    end 
end 

% remove genesis beyond 365 days  sum(sum(kernel(:,:,74)))=5.1735e-04
kernel(:,:,74) = 0;

%normalize the kernel to be 1, so that plot will show correct scale
kernel = kernel/sum(kernel(:));


%% plot the coastline with example density distribution 
long_plot = linspace(250,350,size_lon);
lat_plot = linspace(0,60,size_lat);
time_plot = linspace(0, 365,size_days);
%surface plot
imagesc(long_plot, lat_plot, [kernel(:,:,37)])
hold on;
colorbar;
plot_coast(250:350, 5:50,[.2 .2 .2]);
set(gca,'Ydir','Normal')
grid off

%Consider doing 
%scatter(array_3d_genesis_count(

% [x,y] = meshgrid(long_plot,lat_plot);
% z=squeeze([kernel(:,:,37)]);
% surf(x,y,z);

%% generate number N of genesis
%4-D plot
%sliceomatic([kernel(:,:,37)],long_plot,lat_plot,time_plot);
%sampling    sum(sum(sum(kernel))) = 698.7
pmf = kernel./sum(sum(sum(kernel)));
  %Let N be the desired number of random samples and X, Y, Z as you computed them. Then do:
  N = 50000;
[~,ix] = histc(rand(N,1),cumsum([0;pmf(:)]));
%plot the cummulative density function
x_axes_cumsum = 1:1:size(cumsum([0;pmf(:)]));
line(x_axes_cumsum, cumsum([0;pmf(:)]));
%sample should be random number between two boundaries !!!!!
%sample = [array_3d_genesis_lon(ix),array_3d_genesis_lat(ix),array_3d_genesis_day(ix)];
%within each cell, randmly select a point in the lat, long, time window
sample = zeros(N,3);
for i = 1:N %size(ix)
    %            lon array right boundary  - random number from 0-1, multiply increment   
    sample_lat = array_3d_genesis_lat(ix(i))- rand*0.5;
    sample_lon = array_3d_genesis_lon(ix(i))- rand*0.5;
    sample_day = array_3d_genesis_day(ix(i))- rand*5;
    sample (i, :) = [sample_lat, sample_lon, sample_day];
end

save (['sample_' num2str(N)], 'sample');

scatter3(sample(:,1),sample(:,2),sample(:,3));
hist(sample(:,1));
hist(sample(:,2));
hist(sample(:,3));
%location plot
scatter(sample(:,2),sample(:,1));

