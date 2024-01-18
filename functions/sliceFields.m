function sliceFields(mygrids,meanField,slices,axisrange,n,power,vargin)

%%
mycolormap2 = mycolor('#0000B2','#FFFFFF','#B10000');

if n==1 
    if power ==1
        labelstring1 = '$\langle V_x \rangle (mm/s)$';
        labelstring2 = '$\langle V_y \rangle (mm/s)$';
        labelstring3 = '$\langle V_z \rangle (mm/s)$';
    elseif power ==2
        labelstring1 = '$\sqrt{ \langle V_x^2 \rangle  } (mm/s)$';
        labelstring2 = '$\sqrt{ \langle V_y^2 \rangle  } (mm/s)$';
        labelstring3 = '$\sqrt{ \langle V_z^2 \rangle  } (mm/s)$';
    end
elseif n==2
    if power ==1
        labelstring1 = '$\langle A_x \rangle (mm/s^2)$';
        labelstring2 = '$\langle A_y \rangle (mm/s^2)$';
        labelstring3 = '$\langle A_z \rangle (mm/s^2)$';
    elseif power ==2
        labelstring1 = '$\sqrt{ \langle A_x^2 \rangle } (mm/s^2)$';
        labelstring2 = '$\sqrt{ \langle A_y^2 \rangle } (mm/s^2)$';
        labelstring3 = '$\sqrt{ \langle A_z^2 \rangle } (mm/s^2)$';
    end
end

ifslip =0;
if nargin>6
    % add a particle to plot
    ifslip = vargin(1);
end

if ifslip == 1
    idx = find(sqrt(mygrids.XX.^2+mygrids.YY.^2+mygrids.ZZ.^2)<0.5);
    meanField.x(idx) = NaN;
    meanField.y(idx) = NaN;
    meanField.z(idx) = NaN;
end


%%
figure;
tiledlayout(1,3,'TileSpacing','loose')
nexttile
hf = slice(mygrids.XX,mygrids.YY,mygrids.ZZ,meanField.x,slices.x,slices.y,slices.z);
shading flat;axis equal padded;box
axis(axisrange)
set(gca,FontSize=15)
% title('$Mean\ Field$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
zlabel('$z/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2);
col =colorbar;
if power ==1
%     caxis([-max(abs(xlim(col))),max(abs(xlim(col)))])
end
xlabel(col,labelstring1,'interpreter','latex',FontWeight='bold',FontSize=18)
% xlim(col,[-1e-2 1e-2])


nexttile
hf = slice(mygrids.XX,mygrids.YY,mygrids.ZZ,meanField.y,slices.x,slices.y,slices.z);
shading flat;axis equal padded;box
axis(axisrange)
set(gca,FontSize=15)
% title('$Mean\ Field$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
zlabel('$z/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2);
col =colorbar;
if power ==1
%     caxis([-max(abs(xlim(col))),max(abs(xlim(col)))])
end
xlabel(col,labelstring2,'interpreter','latex',FontWeight='bold',FontSize=18)


nexttile
hf = slice(mygrids.XX,mygrids.YY,mygrids.ZZ,meanField.z,slices.x,slices.y,slices.z);
shading flat;axis equal padded;box
axis(axisrange)
set(gca,FontSize=15)
% title('$Mean\ Field$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
zlabel('$z/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
colormap(mycolormap2);
col =colorbar;
if power ==1
%     caxis([-max(abs(xlim(col))),max(abs(xlim(col)))])
end
xlabel(col,labelstring3,'interpreter','latex',FontWeight='bold',FontSize=18)
