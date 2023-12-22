function quiverFields(XX,YY,ZZ,x,y,z,axisrange)

figure;
quiver3(XX,YY,ZZ,x,y,z,2,LineWidth=1)


xlabel('y');ylabel('x');zlabel('z')
axis(axisrange)
axis equal padded;box
set(gca,FontSize=15)
% title('$Mean\ Field$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$y(g)/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$x/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
zlabel('$z/mm$','interpreter','latex',FontWeight='bold',FontSize=18)
% colormap(mycolormap2);
% col =colorbar;
% xlabel(col,labelstring1,'interpreter','latex',FontWeight='bold',FontSize=18)