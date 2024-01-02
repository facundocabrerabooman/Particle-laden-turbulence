function savefig_custom(fpath,name,x_width,y_width,vargin)

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

set(gcf,'PaperSize',[x_width y_width]); %set the paper size to what you want
%set(gca,'fontsize',14) % ONE FIGURE
%set(gca,'fontsize',17) % TWO FIGURES SIDE BY SIDE

% set(gca,'fontsize',15) 


%%
fname = [fpath filesep name];
savefig(fname)
saveas(gcf,fname,'pdf')
if nargin<5
    saveas(gcf,fname,'png')
else
    format = vargin(1);
    saveas(gcf,[fname '.' format])
end

disp('saved!')

