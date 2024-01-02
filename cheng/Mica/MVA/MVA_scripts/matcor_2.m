function R=matcor(xR,N)

% R=matcor(xR,N)
% estime la matrice de covariance du signal en considérant des fenêtres de N points
% xR est un vecteur
% (c) nmordant dernière modif 27-05-2000

[Mx,Nx]=size(xR);
if Mx>Nx 
   xR=xR.';
end;

Lx=length(xR);
nb=Lx-N+1;	
jkl=1:nb;
nn=(0:N-1)';
jkl=jkl(ones(N,1),:)+ nn(:,ones(size(jkl)));
xR = reshape( xR(jkl), size(jkl) );
R=(xR*xR'+conj(fliplr(flipud(xR)))*fliplr(flipud(xR)).')/2/nb;
