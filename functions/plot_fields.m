function []=plot_fields(folderin,folderout, outputFileName)

load([folderin filesep outputFileName])

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 3D arrows
name = '3D arrows' ;
figure;
q=quiver3(X,Y,Z,mXdt,mYdt,mZdt,'r','LineWidth',2); hold on
[startX,startY,startZ] = meshgrid(5,-25:5:25,-25:5:25);
streamline(X,Y,Z,mXdt,mYdt,mZdt,startX,startY,startZ)

%%%%
%// Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));

%// Get the current colormap
currentColormap = colormap(gca);

%// Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

title(name)
axis equal
q.Marker = '.';
box
xlabel('x')
ylabel('y')
zlabel('z')
savefig_FC([folderout name],8,6,'fig')
savefig_FC([folderout name],8,6,'pdf')

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Get center vertical plane 3D
name = 'center vertical plane3D';
figure;
q=quiver3(X(:,5,:),Y(:,5,:),Z(:,5,:),mXdt(:,5,:),mYdt(:,5,:),mZdt(:,5,:),'r'); hold on
%[startX,startY,startZ] = meshgrid(5,-25:5:25,-25:5:25);
%streamline(X,Y,Z,mXdt,mYdt,mZdt,startX,startY,startZ)

title(name)
axis equal
q.Marker = '.';
box
xlabel('x')
ylabel('y')
zlabel('z')
savefig_FC([folderout name],8,6,'fig')
savefig_FC([folderout name],8,6,'pdf')


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Get center horizontal plane 3D
name = 'center horizontal plane3D';
figure;
q=quiver3(X(:,:,5),Y(:,:,5),Z(:,:,5),mXdt(:,:,5),mYdt(:,:,5),mZdt(:,:,5),'r'); hold on
%[startX,startY,startZ] = meshgrid(5,-25:5:25,-25:5:25);
%streamline(X,Y,Z,mXdt,mYdt,mZdt,startX,startY,startZ)

title(name)
axis equal
q.Marker = '.';
box
xlabel('x')
ylabel('y')
zlabel('z')
savefig_FC([folderout name],8,6,'fig')
savefig_FC([folderout name],8,6,'pdf')

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Get center horizontal plane 3D
close all
name = 'center vertical plane2DY vertvelcolormap';
figure;

q=quiver(Y(:,5,:),Z(:,5,:),mYdt(:,5,:),mZdt(:,5,:),'r','LineWidth',2); hold on

vel = reshape(mZdt(:,5,:),[8,8]);
Ys = reshape(Y(:,5,:),[8,8]);
Zs = reshape(Z(:,5,:),[8,8]);
surf(Ys,Zs,vel)

title('Vertical vel. colormap -- center vertical plane (Y)')
colorbar;caxis([-5.5e-3 5.5e-3]);
axis equal
q.Marker = '.';
xlabel('x')
ylabel('z')
view(0,80)
savefig_FC([folderout name],8,6,'fig')
savefig_FC([folderout name],8,6,'pdf')


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Get center vertical plane 2D Y
name = 'center vertical plane2DY horvelcolormap';
figure;

q=quiver(Y(:,5,:),Z(:,5,:),mYdt(:,5,:),mZdt(:,5,:),'r','LineWidth',2); hold on

vel = reshape(mYdt(:,5,:),[8,8]);
Ys = reshape(Y(:,5,:),[8,8]);
Zs = reshape(Z(:,5,:),[8,8]);
surf(Ys,Zs,real(vel))

title('Hor vel. colormap -- center vertical plane (Y)')
colorbar ;caxis([-5.5e-3 5.5e-3]);
axis equal
q.Marker = '.';
xlabel('x')
ylabel('z')
view(0,80)
savefig_FC([folderout name],8,6,'fig')
savefig_FC([folderout name],8,6,'pdf')



%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Get center vertical plane 2D X
name = 'center vertical plane2DX vertvelcolormap';
figure;

%q=quiver(X(5,:,:),Z(5,:,:),mXdt(5,:,:),mZdt(5,:,:),'r','LineWidth',2); hold on

vel = reshape(mZdt(5,:,:),[8,8]);
Xs = reshape(X(5,:,:),[8,8]);
Zs = reshape(Z(5,:,:),[8,8]);
surf(Xs,Zs,vel)

title('Vertical vel. colormap -- center vertical plane (X fixed)')
colorbar ;caxis([-5.5e-3 5.5e-3]);
axis equal
q.Marker = '.';
xlabel('x')
ylabel('z')
view(0,80)
savefig_FC([folderout name],8,6,'fig')
savefig_FC([folderout name],8,6,'pdf')


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Get center horizontal plane 2D X
name = 'center vertical plane2DX horvelcolormap';
figure;

q=quiver(X(5,:,:),Z(5,:,:),mXdt(5,:,:),mZdt(5,:,:),'r','LineWidth',2); hold on

vel = reshape(mXdt(5,:,:),[8,8]);
Xs = reshape(X(5,:,:),[8,8]);
Zs = reshape(Z(5,:,:),[8,8]);
surf(Xs,Zs,vel)

title('Hor vel. colormap -- center vertical plane (X fixed)')
colorbar ;caxis([-5.5e-3 5.5e-3]);
axis equal
q.Marker = '.';
xlabel('x')
ylabel('z')
view(0,90)
savefig_FC([folderout name],8,6,'fig')
savefig_FC([folderout name],8,6,'pdf')

end