function [mXdt mBdt bins] = track2meanDxDt3DProfile(traj,fieldin,dt,varargin)

% [mXdt mBdt bins] = track2meanDxDt3DProfile(traj,fieldin,dt,varargin)
%
% Calculates the mean profile of (<(d^n X / dt^n)^p> = C*t^(n+p)) + B(noise), using the dt method from a
% Lagrangian trajectories structure. The retrieve quantity mXdt is C^(1/p)
%
% For instance  n = 1 ; p = 1 : C profile of mean velocity
%               n = 1 ; p = 2 : C^1/2 profile of rms velocity
%               n = 2 ; p = 1 : C is the profile of mean acceleration
%               n = 2 ; p = 2 : C^1/2 is the profile of rms acceleration
%               etc.
%
% By default n = 1 and p = 1
%
% Input parameters :
%   traj :      structure of trajectories
%   fieldin :   the field of traj to be differentaited
%   dt :        vector of dt. As a centered difference scheme is used, dt must be a
%               vector of even numbers.
%   varargin{1} : [Nbinsa Nbinsb], number of bins for the profile
%   varargin{2} : order n of the derivative to be estimated (default n = 1)
%   varargin{3} : power p of the averaged quantity (default p = 1 : simple
%                   average)
%   varargin{4} : component of the differentiated quantity
%                ('x','y','z','r','th','z'). If omitted fieldin is used as
%                input coordinnate
%   varargin{5} : system of coordinates cartesian o polar ('cart','pol')
%
%   Output parameters :
%       mXdt : C^1/p estimate of the denoised <(d^n X/dt^n)^p>^1/p
%       mBdt : B^1/p estimate of the noise
%       bins : bins of the profile
%
% MB 28/04/2020
%


o = 1;
p = 1;
power = 1;
system = 'cart';
coord = fieldin;
Nbina = 16;
Nbinb = 17;
Nbinc = 18;

if nargin > 3
    Nbina=varargin{1}(1);
    Nbinb=varargin{1}(2);
    Nbinc=varargin{1}(3);
end
if nargin > 4
    o = varargin{2};
end
if nargin > 5
    power = varargin{3};
end
if nargin > 6
    coord = varargin{4};
end
if nargin > 7
    system = varargin{5};
end


traj = addStructFun(traj,'X','N',@(X)(numel(X)));
NN = [traj.N];
Iinc = find(NN>2*max(dt));
traj = traj(Iinc);
%%

trajinc = addStructFun(traj,'X','xr',@(X)(X(1+o*max(dt)/2:end-o*max(dt)/2)));
trajinc = addStructFun(trajinc,'Y','yr',@(X)(X(1+o*max(dt)/2:end-o*max(dt)/2)));
trajinc = addStructFun(trajinc,'Z','zr',@(X)(X(1+o*max(dt)/2:end-o*max(dt)/2)));


switch coord
    case{'x'}
        trajinc = addStructFun(trajinc,'X','dX',@(X)(finiteCenteredDiffn(X,dt,o,p)));
    case{'y'}
        trajinc = addStructFun(trajinc,'Y','dX',@(X)(finiteCenteredDiffn(X,dt,o,p)));
    case{'z'}
        trajinc = addStructFun(trajinc,'Z','dX',@(X)(finiteCenteredDiffn(X,dt,o,p)));
    case{'r'}
        trajinc = addStructFun(trajinc,'X','dx',@(X)(finiteCenteredDiffn(X,dt,o,p)));
        trajinc = addStructFun(trajinc,'Y','dy',@(X)(finiteCenteredDiffn(X,dt,o,p)));
        for k = 1:numel(trajinc)
            [th,r,z] = cart2pol(trajinc(k).xr,trajinc(k).yr,trajinc(k).zr);
            trajinc(k).r = r;
            trajinc(k).th = th;
            trajinc(k).dX = trajinc(k).dx .* cos(trajinc(k).th) +  trajinc(k).dy .* sin(trajinc(k).th);
        end
    case{'th'}
        trajinc = addStructFun(trajinc,'X','dx',@(X)(finiteCenteredDiffn(X,dt,o,p)));
        trajinc = addStructFun(trajinc,'Y','dy',@(X)(finiteCenteredDiffn(X,dt,o,p)));
        for k = 1:numel(trajinc)
            [th,r,z] = cart2pol(trajinc(k).xr,trajinc(k).yr,trajinc(k).zr);
            trajinc(k).r = r;
            trajinc(k).th = th;
            trajinc(k).dX = -trajinc(k).dx .* sin(trajinc(k).th) +  trajinc(k).dy .* cos(trajinc(k).th);
        end
    otherwise
        trajinc = addStructFun(trajinc,fieldin,'dX',@(X)(finiteCenteredDiffn(X,dt,o,p)));
        
end

dX = vertcat(trajinc.dX);

a = vertcat(trajinc.xr);
b = vertcat(trajinc.yr);
c = vertcat(trajinc.zr);

switch system
    case{'polar'}
        [a,b,c] = cart2pol(a,b,c);
end
%%

edgesa = linspace(min(a)-1e-3*abs(min(a)),max(a)+1e-3*abs(max(a)),Nbina+1);
edgesb = linspace(min(b)-1e-3*abs(min(b)),max(b)+1e-3*abs(max(b)),Nbinb+1);
edgesc = linspace(min(c)-1e-3*abs(min(c)),max(c)+1e-3*abs(max(c)),Nbinc+1);
[mdX edges mid] = arrayfun(@(I)(histcn([a b c],edgesa,edgesb,edgesc,'AccumData',dX(:,I),'Fun',@(X)(mean(X.^power)))),1:numel(dt),'UniformOutput',false);
bins = mid{1};

pp = polyfithypercube(dt.^(power+o),mdX,1);
mXdt = cellfun(@(C)(C(1)),pp).^(1/(power));
mBdt = cellfun(@(C)(C(2)),pp).^(1/(power));
