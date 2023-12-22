function grille=calcParamGrille(x,U)

%
% grille=paramGrille(x,U)
% 
% calcule les parametres de la turbulence de la grille a une distance x, a
% partir des mesures fil chaud dans la these de nicolas Mazellier
%
% MB 2007
%
M=7.5e-2;
%M=.3;
x0=2*M;
%x0=0;
nu=1.5e-5;

grille.M=M;
grille.x0=x0;
grille.nu=nu;
grille.U=U;
grille.x=x;

grille.tau=sqrt(0.073*((x-x0)/M).^-1.17);
grille.L=0.041*(x-x0).^.40;
grille.eta=9.6e-4*U.^(-3/4)*(x-x0).^0.61*(7.5e-2/M)^0.61;

grille.Rlambda=sqrt(15.*grille.tau.*U.*grille.L/nu);
grille.lambda=15*grille.L./grille.Rlambda;

grille.epsilon=nu^3./grille.eta.^4;

grille.tau_eta=sqrt(nu./grille.epsilon);
grille.Te=grille.L./grille.tau./U;