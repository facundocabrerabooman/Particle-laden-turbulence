function trajinc = dtVel(traj,fieldin,fieldout,dt,power)
% 
% trajinc = dtVel(traj,fieldin,fieldout,Ninc,power)
%
% estimates the local 1st order centered increment with dt = Ninc of (traj.(fieldin)^power) at each trajectory position 
%

traj = addStructFun(traj,fieldin,'N',@(X)(numel(X)));

NN = [traj.N];
Iinc = find(NN>max(dt));

trajinc = addStructFun(traj(Iinc),fieldin,'dx',@(X)(finiteCenteredDiffn(X,dt,power)));

%%
%F = arrayfun(@(X)((X.dx{1}(1:X.N-5)+X.dx{2}(1:X.N-5)/2+X.dx{3}(1:X.N-5)/3+X.dx{4}(1:X.N-5)/4+X.dx{5}(1:X.N-5)/5)/5),traj,'UniformOutput','false');

%%
for k=1:numel(trajinc)
    dx = trajinc(k).dx./dt;
    trajinc(k).(fieldout) = mean(dx,2);    
end

%%