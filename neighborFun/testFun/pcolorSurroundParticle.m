function pcolorSurroundParticle(neighbor_terminalState,field)

Rmin = neighbor_terminalState(1).Rmin;
Rmax = neighbor_terminalState(1).Rmax;
dR   = neighbor_terminalState(1).dR;

n = (Rmax-Rmin)/dR;

% Define polar coordinates (theta and radius)
nc = 4;
theta = linspace(0, 360, 100*nc);
radius = linspace(Rmin, Rmax, n);  % Adjust the radius range as needed

for i = 1:numel(radius)
    if strcmp(field,'V')
        frontx(1:numel(theta)/nc,i) = neighbor_terminalState(i).VrelFrontx;
        fronty(1:numel(theta)/nc,i) = neighbor_terminalState(i).VrelFronty;
        frontz(1:numel(theta)/nc,i) = neighbor_terminalState(i).VrelFrontz;
        front (1:numel(theta)/nc,i) = neighbor_terminalState(i).VrelFront;
        backx(1:numel(theta)/nc,i)  = neighbor_terminalState(i).VrelBackx;
        backy(1:numel(theta)/nc,i)  = neighbor_terminalState(i).VrelBacky;
        backz(1:numel(theta)/nc,i)  = neighbor_terminalState(i).VrelBackz;
        back (1:numel(theta)/nc,i)  = neighbor_terminalState(i).VrelBack;
    elseif strcmp(field,'A')
        frontx(1:numel(theta)/nc,i) = neighbor_terminalState(i).ArelFrontx;
        fronty(1:numel(theta)/nc,i) = neighbor_terminalState(i).ArelFronty;
        frontz(1:numel(theta)/nc,i) = neighbor_terminalState(i).ArelFrontz;
        front (1:numel(theta)/nc,i) = neighbor_terminalState(i).ArelFront;
        backx(1:numel(theta)/nc,i)  = neighbor_terminalState(i).ArelBackx;
        backy(1:numel(theta)/nc,i)  = neighbor_terminalState(i).ArelBacky;
        backz(1:numel(theta)/nc,i)  = neighbor_terminalState(i).ArelBackz;
        back (1:numel(theta)/nc,i)  = neighbor_terminalState(i).ArelBack;
    end
end

figure;
polarPcolor(radius,theta, data);
colormap jet;
colorbar;
