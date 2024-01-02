function [myMeanSlice,sliceData] = sliceSlipVelo(grids,meanFields,nSlice,rp,Rmax,nSpacing,n,power,vargin)

if nargin>8
    viewSlices = vargin(1);
end

fieldname = fieldnames(meanFields);
for i = 1:numel(fieldname)
    [myMeanSlice.(fieldname{i}),sliceData.(fieldname{i})]  = meanSlice(grids,meanFields.(fieldname{i}),nSlice,Rmax,nSpacing);
end

%% here we switch x and y axis
% xData = sliceData.x(1).xData(1,:);
% yData = sliceData.x(1).yData(:,1);
% 
% % set index inside the particles to NaN
% d = sqrt(xData.^2+yData.^2);
% for i = 1:numel(fieldname)
%     myMeanSlice.(fieldname{i})(d<rp) = NaN;    
% end

xData = sliceData.x(1).xData;
yData = sliceData.x(1).yData;
zData = zeros(size(xData));

d = sqrt(xData.^2+yData.^2);
for i = 1:numel(fieldname)
    myMeanSlice.(fieldname{i})(d<rp) = NaN;    
end

% Transpose the matrix to swith the x and y axis
for i = 1:numel(fieldname)
    cData.(fieldname{i}) = (myMeanSlice.(fieldname{i}))';
end

%%
mycolormap2 = mycolor('#0000B2','#FFFFFF','#B10000');

if n==1 
    if power ==1
        labelstring1 = '$\langle V^{slip}_x \rangle (mm/s)$';
        labelstring2 = '$\langle V^{slip}_y \rangle (mm/s)$';
        labelstring3 = '$\langle V^{slip}_z \rangle (mm/s)$';
    elseif power ==2
        labelstring1 = '$\sqrt{ \langle (V^{slip}_x)^2 \rangle  } (mm/s)$';
        labelstring2 = '$\sqrt{ \langle (V^{slip}_y)^2 \rangle  } (mm/s)$';
        labelstring3 = '$\sqrt{ \langle (V^{slip}_z)^2 \rangle  } (mm/s)$';
    end
elseif n==2
    if power ==1
        labelstring1 = '$\langle A^{slip}_x \rangle (mm/s^2)$';
        labelstring2 = '$\langle A^{slip}_y \rangle (mm/s^2)$';
        labelstring3 = '$\langle A^{slip}_z \rangle (mm/s^2)$';
    elseif power ==2
        labelstring1 = '$\sqrt{ \langle (A^{slip}_x)^2 \rangle } (mm/s^2)$';
        labelstring2 = '$\sqrt{ \langle (A^{slip}_x)^2 \rangle } (mm/s^2)$';
        labelstring3 = '$\sqrt{ \langle (A^{slip}_x)^2 \rangle } (mm/s^2)$';
    end
end


%% plot
figure;
tiledlayout(1,3,'TileSpacing','loose')
nexttile
% hf = imagesc('XData',xData,'YData',yData,'CData',cData.(fieldname{1}));
surf(xData, yData,zData,cData.(fieldname{1}))
view(0,90)
shading interp;axis equal tight;box
set(gca,FontSize=15)
xlabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2);
col =colorbar;
if power ==1
    caxis([-max(abs(xlim(col))),max(abs(xlim(col)))])
end
xlabel(col,labelstring1,'interpreter','latex',FontWeight='bold',FontSize=18)
% xlim(col,[-1e-2 1e-2])


nexttile
% hf = imagesc('XData',xData,'YData',yData,'CData',cData.(fieldname{2}));
surf(xData, yData,zData,cData.(fieldname{2}))
view(0,90)
shading interp;axis equal tight;box
set(gca,FontSize=15)
xlabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2);
col =colorbar;
if power ==1
    caxis([-max(abs(xlim(col))),max(abs(xlim(col)))])
end
xlabel(col,labelstring2,'interpreter','latex',FontWeight='bold',FontSize=18)


nexttile
% hf = imagesc('XData',xData,'YData',yData,'CData',cData.(fieldname{3}));
surf(xData, yData,zData,cData.(fieldname{3}))
view(0,90)
shading interp;axis equal tight;box
set(gca,FontSize=15)
xlabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2);
col =colorbar;
if power ==1
    caxis([-max(abs(xlim(col))),max(abs(xlim(col)))])
end
xlabel(col,labelstring3,'interpreter','latex',FontWeight='bold',FontSize=18)

%%
if viewSlices
    figure
    for i  = 1:numel(sliceData.x)
        surf(sliceData.x(i).xData,sliceData.x(i).yData,sliceData.x(i).zData,sliceData.x(i).cData)
        hold on
    end
    xlabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
    ylabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
    zlabel('$z/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
end

