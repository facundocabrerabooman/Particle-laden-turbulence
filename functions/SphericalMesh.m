function [r, theta, rgrid, tgrid, X, Y] = SphericalMesh(radius, nbinr, nbint, theta)
% SPHERICALMESH generates a spherical mesh grid in Cartesian coordinates.
%   [r, theta, rgrid, tgrid, X, Y] = SPHERICALMESH(radius, nbinr, nbint, theta)
%   creates a spherical mesh grid with specified radial and azimuthal coordinates.
%
%   Input:
%       radius - Vector [Rmin, Rmax] specifying the radial bounds.
%       nbinr  - Number of bins in the radial direction.
%       nbint  - Number of bins in the azimuthal direction.
%       theta  - Optional parameter specifying the azimuthal angle range.
%
%   Output:
%       r      - Radial coordinates.
%       theta  - Azimuthal coordinates.
%       rgrid  - Mesh grid of radial coordinates.
%       tgrid  - Mesh grid of azimuthal coordinates.
%       X      - Cartesian X coordinates.
%       Y      - Cartesian Y coordinates.
%

    % Default values for nbinr and nbint if not provided
    if nargin < 2
        nbinr = 20;
        nbint = 20;
    end
    
    % Default azimuthal angle range if not provided
    if nargin < 4
        thetaMin = 0;
        thetaMax = pi;
    else
        thetaMin = theta(1);
        thetaMax = theta(2);
    end
    
    % Extract radial bounds
    Rmin = radius(1);
    Rmax = radius(2);

    % Define spherical grid parameters
    r = linspace(Rmin, Rmax, nbinr); % radial coordinates
    theta = linspace(thetaMin, thetaMax, nbint); % azimuthal coordinates

    % Create the mesh grid in spherical coordinates
    [rgrid, tgrid] = meshgrid(r, theta);

    % Convert spherical coordinates to Cartesian coordinates
    X = rgrid .* cos(tgrid);
    Y = rgrid .* sin(tgrid);

    % Uncomment the following lines to visualize the 2D mesh with grid lines
%     figure;
%     plot(X, Y, '.');
%     hold on;
%     for i = 1:size(X, 2)
%         plot(X(:, i), Y(:, i), '-k');
%     end
%     for i = 1:size(X, 1)
%         plot(X(i, :), Y(i, :), '-k');
%     end
%     hold off;
%     title('Spherical Mesh Grid (2D) with Grid Lines');
%     xlabel('X');
%     ylabel('Y');
%     axis equal;
%     grid on;
end

