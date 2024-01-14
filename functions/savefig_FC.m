function []=savefig_FC(name,x_width,y_width, format,saveflag)

set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperPosition', [0.1 0.27 x_width y_width]); %
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

set(gcf,'PaperSize',[x_width y_width]); %set the paper size to what you want  
set(gca,'fontsize',20) % one figure, expif ddt exp

%set(gca,'fontsize',20) % RSI DDT multiphase

%%%%
if ~exist('saveflag')
 if format == 'pdf'
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(gcf,name,'-painters','-dpdf','-r0')
 saveas(gcf,[name '.' format]) 
 else 
     saveas(gcf,[name '.' format]) 
 end
end 
end

