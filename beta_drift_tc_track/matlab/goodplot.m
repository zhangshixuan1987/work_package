function goodplot(papersizeW,papersizeH,margin,fontsize)
% function which produces a nice-looking plot
% and sets up the page for nice printing
% Adi Nugraha
% PNNL 05/2015


%set(gca, 'Position',[margin margin 1 1])

cvt=2.54;

if nargin == 0
margin = 0.01;
fontsize = 10;
elseif nargin == 2
margin = 0.01;
fontsize = 10;
elseif nargin == 3
fontsize = 10;
end

papersizeW=papersizeW/cvt;
papersizeH=papersizeH/cvt;
margin=margin/cvt;

set(get(gca,'xlabel'),'FontSize', fontsize, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fontsize, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', 13, 'FontWeight', 'Bold');
%set(get(gca,'xlabel'),'FontSize', fontsize);
%set(get(gca,'ylabel'),'FontSize', fontsize);
%set(get(gca,'title'),'FontSize', fontsize);
%set(gca,'XMinorTick','on','YMinorTick','on');
%set(gca,'ticklength',2*get(gca,'ticklength'))
%box off; 
%axis square;
set(gca,'LineWidth',1);
set(gca,'FontSize',fontsize);
set(gca,'TickDir','out');
%set(gca,'FontWeight','Bold');
set(gcf,'color','w');

%# set size of figure's "drawing" area on screen
set(gcf, 'Units','inches', 'Position',[0.5/cvt 0.5/cvt papersizeW-2*margin papersizeH-2*margin])

set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize', [papersizeW papersizeH]);
set(gcf,'PaperPosition',[margin margin papersizeW-2*margin papersizeH-2*margin]);
%set(gcf,'PaperPositionMode','Manual');

%# WYSIWYG mode: you need to adjust your screen's DPI (*)
%set(gcf, 'PaperPositionMode','auto')
