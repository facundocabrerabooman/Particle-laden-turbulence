function params = expParam(Fs,mp,dp,rhof,nu,geff,taup,Vs,epsilon,sigmau)
% compute turbulence quantities
% [ly,lx,Gamma,beta,taup0,Vg] = turb_param(mp,dp,rhof,nu,px,py,geff)
%
% INPUT:
% rhop: particle mass       (kg)
% dp: particle diameter     (m)
% rhof: fluid density       (kg/m3)
% nu: fluid viscosity       (m2/s)
% px: pixel size of image
% py: pixel size of image
% geff: ratio of the effective gravity to 9.8m/s2 
% taup: estunated from (V(t)-V0)/(Vs-V0)
% Vs: measured particle's terminal falling velocity (m/s)
% epsilon: disspation rate
% sigmau: rms value of the fluctuation velocity of fluid 
%
% OUTPUT:
% ly: domain length along gravity (vertical)                (m)
% lx: domain length in horizontal direction                 (m)
% Gamma: density ratio
% beta: 3*rho_f/(rho_f+2*rho_p)
% taup0: particle response time with linear approximation   (s)
% Vg: terminal settling velocity with linear approximation  (s)
%
% reference turbulence parameters
% mp = 0.0048e-3;%kg
% dp = 1.0094e-3;%m
% rhof = 1000; %kg/m3
% nv = 8.532e-7;%m2/s, @ 27 celsius
% geff = 1; %ratio to 9.8m/s2
% turb = turb_param(0.0048e-3,1.0094e-3,1000,8.532e-7,1)

%% measured values
params.Fs = Fs;
params.mp = mp;
params.dp = dp;
params.rhop = params.mp/(pi*4/3*(params.dp/2)^3);
params.rhof = rhof;
params.geff = geff*9.8;
params.nu = nu;
params.taup = taup;
params.Vs = Vs; 


% % calibration plate type: 058-5
% mire_z = 5e-3; % m
% mire_pixel = 85;
% % fps = 3000;
% delta_pixel = 10; % 
% 
% % py = 1280;
% % px = 512;
% 
% % exp param
% sensor_resol = mire_z/mire_pixel; % m/px
% ly = py*sensor_resol; % m
% lx = px*sensor_resol; % m

%% real parameters from measurements
% density ratio
params.Gamma = params.rhop/params.rhof; 
params.beta = 3/(1+2*params.Gamma); % 3*rho_f/(rho_f+2*rho_p)

% measured particle Reynolds number at terminal state
params.Reps_real = Vs*dp/nu;

% CD from measured taup
params.CD_real = 4*params.rhop*params.dp^2/(3*params.rhof*params.taup*params.nu)*params.Rep^(-1);

%% fluid parameters:
% kolmogrov scale
params.eta = (nu^3/epsilon);

% Taylor scale
params.lambda = sqrt(15*nu/epsilon)*sigmau;

% kolmogrov time scale
params.tau_eta = sqrt(nu/epsilon);

% Taylor scale Reynolds number
params.Re_lambda = sigmau*params.lambda/nu;

% turbulent instensity
params.Ti = sigmau/Vs;


%% parameters from approximation 
% particle response time with linear approximation
params.taup_LinearApproxi = params.dp^2/(12*params.nu*params.beta);

% teminal velocity from linear approximation
params.Vs_LinearApproxi = (params.Gamma-1)/(params.Gamma+0.5)*params.taup0*params.geff; 

% teminal particle Reynolds number from linear approximation
params.Reps_LinearApproxi = params.Vs*params.dp/params.nu;

% buoyancy velocity and Galileo number
params.Ug = sqrt(abs(params.Gamma-1)*params.geff*params.dp);
params.Galileo = params.Ug*params.dp/params.nu;

% Facu's arxiv: emprical relation approximation
params.CD_GalileoApproxi = 24/params.Rep*(1+0.150*params.Rep^0.681)+0.407/(1+8710/params.Rep);
params.Reps_GalileoApproxi = params.Galileo^2*(22.5+params.Galileo^1.364)/(0.0258*params.Galileo^2.6973+2.81*params.Galileo^2.0306+18*params.Galileo^1.364+405);





