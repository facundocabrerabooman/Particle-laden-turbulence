function vel=traite_serie(serie,ideb,iend,thresh,varargin)

%dirname='/Users/mbourgoi/Documents/Recherche/Lagrangien/Analyse_donnees/051027';
%sname=sprintf('%s/serie%i_vel_N7_K13.mat',dirname,serie);
%fname=sprintf('serie%i_80k.%04i',serie,ideb);
%thresh=5.e-5;
%dirname='d:\qureshi\Lagrangien\051221';
dirname='.';
sname=sprintf('%s/serie%i_vel_N7_K13.mat',dirname,serie);
fname=sprintf('serie%i.%04i',serie,ideb);

disp(sprintf('Fichier %s (%i/%i)',fname,ideb,iend-ideb+1));

if nargin>4
    vel=varargin{1};
    [vel]=velocities(fname,thresh,vel);
else
    [vel]=velocities(fname,thresh);
end
eval(['save ' sname ' vel;']);


for j=ideb+1:iend
   %fname=sprintf('serie%i_80k.%04i',serie,j);
    fname=sprintf('serie%i.%04i',serie,j);
   
   disp(sprintf('Fichier %s (%i/%i)',fname,j,iend-ideb+1));
   [vel]=velocities(fname,thresh,vel);
   if   ((mod(j,10)==0)||(j==iend))
         eval(['save ' sname ' vel;']);
   end
end
