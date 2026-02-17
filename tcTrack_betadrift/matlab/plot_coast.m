function plot_coast(x,y,color)
% plot the world coast boundary
% Xiuquan Wan
% 2004/9/12
%
if nargin<3 color=0; end;
latmin = min(y);
latmax = max(y);
lonmin = min(x);
lonmax = max(x);

load coast

lat2=[nan;lat;nan; lat];
long2=[nan;long; nan; long+360];

plot(long2,lat2,'k','linewidth',1)
axis([lonmin lonmax round(latmin) round(latmax)]);
grid on
hold on

if color~=0;
   a=long2;
   k=[find(isnan(a))];
   for i=1:length(k)-1,
   x=a([k(i)+1:(k(i+1)-1) k(i)+1]);
   y=lat2([k(i)+1:(k(i+1)-1) k(i)+1]);
%   patch(x,y,[0.7 0.4 0.3]);
   patch(x,y,color)
   end;
end
