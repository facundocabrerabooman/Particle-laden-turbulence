function savefig_custom(name,x_width,y_width, format,saveflag)

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

set(gcf,'PaperSize',[x_width y_width]); %set the paper size to what you want
%set(gca,'fontsize',14) % ONE FIGURE
%set(gca,'fontsize',17) % TWO FIGURES SIDE BY SIDE

set(gca,'fontsize',15) 


%%%%
if ~exist('saveflag')
    if format == 'pdf'
        saveas(gcf,[name '.' format])
    else
        saveas(gcf,[name '.' format])
    end
end
end

