% Program to load the reanalysis data
% load netcdf that Karthik provided. It is for 
%Katrina, 2005
%Ike, 2008
%Sandy, 2012 

clear all;
clc;
% reanalysis data for katrina, ike, and sandy
nc_f = netcdf.open('kis_era.nc','NOWRITE');

%[1]display header of the netcdf
ncdisp('kis_era.nc')

% save to array
kis.lon = double(netcdf.getVar(nc_f, 0, [0],[201]));
kis.lat = double(netcdf.getVar(nc_f, 1, [0],[121]));
kis.level = double(netcdf.getVar(nc_f, 2, [0],[3]));

kis.time = zeros(368,1);  
for t = 1:368
    if t<=124   %katrina
        base= date2doy(datenum('8/1/2005'));%using 2005 as base line
        t_start = 1;
    elseif t<=244 && t>=125   %ike
        base= date2doy(datenum('9/1/2008'))+1095;%using 2005 as base line
        t_start = 125;
    else  % sandy
        base= date2doy(datenum('10/1/2012'))+2556;%using 2005 as base line
        t_start = 245;
    end
    
    kis.time(t) = base+(t-t_start)*6.0/24.0;
end
    
kis.u200 = squeeze(netcdf.getVar(nc_f, 4, [0,0,0,0],[201,121,1,368]));
kis.u600 = squeeze(netcdf.getVar(nc_f, 4, [0,0,1,0],[201,121,1,368]));
kis.u850 = squeeze(netcdf.getVar(nc_f, 4, [0,0,2,0],[201,121,1,368]));

kis.v200 = squeeze(netcdf.getVar(nc_f, 5, [0,0,0,0],[201,121,1,368]));
kis.v600 = squeeze(netcdf.getVar(nc_f, 5, [0,0,1,0],[201,121,1,368]));
kis.v850 = squeeze(netcdf.getVar(nc_f, 5, [0,0,2,0],[201,121,1,368]));

% remove falt values
kis.u200(find(kis.u200==-32767)) = NaN;
kis.u600(find(kis.u600==-32767)) = NaN;
kis.u850(find(kis.u850==-32767)) = NaN;
kis.v200(find(kis.v200==-32767)) = NaN;
kis.v600(find(kis.v600==-32767)) = NaN;
kis.v850(find(kis.v850==-32767)) = NaN;

%apply scale factor and offset
kis.u200 = double(kis.u200) * 0.0022186+17.6861;
kis.u600 = double(kis.u600) * 0.0022186+17.6861;
kis.u850 = double(kis.u850) * 0.0022186+17.6861;
kis.v200 = double(kis.v200) * 0.0023972-1.962;
kis.v600 = double(kis.v600) * 0.0023972-1.962;
kis.v850 = double(kis.v850) * 0.0023972-1.962;

%% prepare lat lon time array for interpolation
[ kis.lat_3d,kis.lon_3d, kis.time_3d] = meshgrid(kis.lat, kis.lon, kis.time);
%% plot to check the distribution of u,v
pcolor(kis.u850(:,:,122));

%% plot example vorticity %skip the boundaries
kis.vort_200= NaN(201,121,368);
kis.vort_600= NaN(201,121,368); 
kis.vort_850= NaN(201,121,368);
for i = 2:200
    for j = 2:120
        for time = 1:368
            kis.vort_200(i,j,time) = (kis.v200(i+1,j,time) - kis.v200(i-1,j,time))/...
                (kis.lon_3d(i+1,j,time) - kis.lon_3d(i-1,j,time)) - ...
                (kis.u200(i,j+1,time) - kis.u200(i,j-1,time))/...
                (kis.lat_3d(i,j+1,time) - kis.lat_3d(i,j-1,time));
                
            kis.vort_600(i,j,time) = (kis.v600(i+1,j,time) - kis.v600(i-1,j,time))/...
                (kis.lon_3d(i+1,j,time) - kis.lon_3d(i-1,j,time)) - ...
                (kis.u600(i,j+1,time) - kis.u600(i,j-1,time))/...
                (kis.lat_3d(i,j+1,time) - kis.lat_3d(i,j-1,time));
            
            kis.vort_850(i,j,time) = (kis.v850(i+1,j,time) - kis.v850(i-1,j,time))/...
                (kis.lon_3d(i+1,j,time) - kis.lon_3d(i-1,j,time)) - ...
                (kis.u850(i,j+1,time) - kis.u850(i,j-1,time))/...
                (kis.lat_3d(i,j+1,time) - kis.lat_3d(i,j-1,time));
        end
    end
end

% plot the vorticity
h = pcolor(kis.lon, kis.lat, kis.vort_850(:,:,110)');
hold on;
set(h,'edgecolor','none')
colorbar;
title('vorticity 850');
plot_coast(260:350, 5:50);
set(gca,'Ydir','Normal')
grid off
print(gcf, '-dpng', sprintf('-r%d',120),'vorticity_850_before');
close;

%% filter the hurricane influenced wind from reanalysis data
load('array_record');
nc_f = netcdf.open('attracks.nc','NOWRITE');
year = netcdf.getVar(nc_f, 0, 0,1753);
month = netcdf.getVar(nc_f, 3, [0,0],[191,1753]);
day = netcdf.getVar(nc_f, 4, [0,0],[191,1753]);
hour = netcdf.getVar(nc_f, 5, [0,0],[191,1753]);

lat_hist = netcdf.getVar(nc_f, 6, [0,0],[191,1753]);
lon_hist = netcdf.getVar(nc_f, 7, [0,0],[191,1753]);

%covert zero to nan
lat_hist(lat_hist == 0) = NaN;
lon_hist(lon_hist < 200) = NaN;
N = 3;
sample_index = [1606,1661,1739];
sample_lon = squeeze(lon_hist(1,sample_index));
sample_lat = squeeze(lat_hist(1,sample_index));
sample_time(1) = date2doy(datenum('8/1/2005'))+day(1,1606)-1+ hour(1,1606)/24.0;
sample_time(2) = date2doy(datenum('9/1/2008'))+1095+day(1,1661)-1+ hour(1,1661)/24.0;
sample_time(3) = date2doy(datenum('10/1/2012'))+2556+day(1,1739)-1+ hour(1,1739)/24.0;

%radius at 1000 mb is 800, radius at 100 mb is 500. Linearly interpolate:
radius = [533, 667, 750];

%define the corresponding track and time steps
kis.track_index = zeros(368,1);
kis.step_index = zeros(368,1);
kis.track_index(92:122) = 1606;
kis.track_index(126:183) = 1661;
kis.track_index(328:367) = 1739;
kis.track = [1,2,3];
kis.step_index(92:122) = 1:31;
kis.step_index(126:183) = 1:58;
kis.step_index(328:367) = 1:40;

%each time step
for time = 1:368
    if kis.step_index(time) ~= 0
        %identify which track this is
        if t<=124   %katrina
            track = 1;
        elseif t<=244 && t>=125   %ike
            track = 2;
        else  % sandy
            track = 3;
        end
        %find center
        lat = lat_hist(kis.step_index(time),kis.track_index(time));
        lon = lon_hist(kis.step_index(time),kis.track_index(time));
        %selet within radius
        for i = 1:201
            for j = 1:121
                lon_1 = kis.lon(i);
                lat_1 = kis.lat(j);
                %
                radius_degree(1) = radius(1)*1000*degree_per_meter( lat, 'lat' );
                radius_degree(2) = radius(2)*1000*degree_per_meter( lat, 'lat' );
                radius_degree(3) = radius(3)*1000*degree_per_meter( lat, 'lat' );
                %200 mb
                if ((lat-lat_1)^2 +(lon-lon_1)^2 )^0.5 <radius_degree(1)
                    %replace to NaN
                    kis.u200(i,j,time) = NaN;
                    kis.v200(i,j,time) = NaN;
                end
                %600 mb
                if ((lat-lat_1)^2 +(lon-lon_1)^2 )^0.5 <radius_degree(2)
                    %replace to NaN
                    kis.u600(i,j,time) = NaN;
                    kis.v600(i,j,time) = NaN;
                end 
                %850 mb
                if ((lat-lat_1)^2 +(lon-lon_1)^2 )^0.5 <radius_degree(3)
                    %replace to NaN
                    kis.u850(i,j,time) = NaN;
                    kis.v850(i,j,time) = NaN;
                end 
            end
        end
        %200 mb
        %interpolate linearly using griddata %count nans
        Inan=find(isnan(kis.u200 )==1 & kis.time_3d ==kis.time(time));
        Igood=find(~isnan(kis.u200)==1 & kis.time_3d ==kis.time(time));
        if size(Igood,1)>0
            kis.u200(Inan)= griddata(kis.lat_3d(Igood),kis.lon_3d(Igood),...
                kis.u200(Igood), kis.lat_3d(Inan),  kis.lon_3d(Inan),'linear'); 
            kis.v200(Inan)= griddata(kis.lat_3d(Igood),kis.lon_3d(Igood),...
                kis.v200(Igood), kis.lat_3d(Inan),  kis.lon_3d(Inan),'linear'); 
        end 
        %600 mb
        %interpolate linearly using griddata %count nans
        Inan=find(isnan(kis.u600 )==1 & kis.time_3d ==kis.time(time));
        Igood=find(~isnan(kis.u600)==1 & kis.time_3d ==kis.time(time));
        if size(Igood,1)>0
            kis.u600(Inan)= griddata(kis.lat_3d(Igood),kis.lon_3d(Igood),...
                kis.u600(Igood), kis.lat_3d(Inan),  kis.lon_3d(Inan),'linear'); 
            kis.v600(Inan)= griddata(kis.lat_3d(Igood),kis.lon_3d(Igood),...
                kis.v600(Igood), kis.lat_3d(Inan),  kis.lon_3d(Inan),'linear'); 
        end 
        %850 mb
        %interpolate linearly using griddata %count nans
        Inan=find(isnan(kis.u850 )==1 & kis.time_3d ==kis.time(time));
        Igood=find(~isnan(kis.u850)==1 & kis.time_3d ==kis.time(time));
        if size(Igood,1)>0
            kis.u850(Inan)= griddata(kis.lat_3d(Igood),kis.lon_3d(Igood),...
                kis.u850(Igood), kis.lat_3d(Inan),  kis.lon_3d(Inan),'linear'); 
            kis.v850(Inan)= griddata(kis.lat_3d(Igood),kis.lon_3d(Igood),...
                kis.v850(Igood), kis.lat_3d(Inan),  kis.lon_3d(Inan),'linear'); 
        end 
    end
end

%% calculate vorticity again and plot
for i = 2:200
    for j = 2:120
        for time = 1:368
            kis.vort_200(i,j,time) = (kis.v200(i+1,j,time) - kis.v200(i-1,j,time))/...
                (kis.lon_3d(i+1,j,time) - kis.lon_3d(i-1,j,time)) - ...
                (kis.u200(i,j+1,time) - kis.u200(i,j-1,time))/...
                (kis.lat_3d(i,j+1,time) - kis.lat_3d(i,j-1,time));
                
            kis.vort_600(i,j,time) = (kis.v600(i+1,j,time) - kis.v600(i-1,j,time))/...
                (kis.lon_3d(i+1,j,time) - kis.lon_3d(i-1,j,time)) - ...
                (kis.u600(i,j+1,time) - kis.u600(i,j-1,time))/...
                (kis.lat_3d(i,j+1,time) - kis.lat_3d(i,j-1,time));
            
            kis.vort_850(i,j,time) = (kis.v850(i+1,j,time) - kis.v850(i-1,j,time))/...
                (kis.lon_3d(i+1,j,time) - kis.lon_3d(i-1,j,time)) - ...
                (kis.u850(i,j+1,time) - kis.u850(i,j-1,time))/...
                (kis.lat_3d(i,j+1,time) - kis.lat_3d(i,j-1,time));
        end
    end
end
% plot the vorticity
h = pcolor(kis.lon, kis.lat, kis.vort_850(:,:,110)');
hold on;
set(h,'edgecolor','none')
colorbar;
title('vorticity 850');
plot_coast(260:350, 5:50);
set(gca,'Ydir','Normal')
grid off
print(gcf, '-dpng', sprintf('-r%d',120),'vorticity_850_after');
close;
%% save data
save('kis.mat', 'kis');





