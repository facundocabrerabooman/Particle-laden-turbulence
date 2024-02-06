function [rho, theta, z] = cartToCylind(X1, X2, RM)
% CARTTOCYLIND converts Cartesian coordinates to cylindrical coordinates.
%   [rho, theta, z] = CARTTOCYLIND(X1, X2, RM) calculates cylindrical
%   coordinates based on input positions X1 and X2, and an optional rotation
%   matrix RM.
%
%   Input:
%       X1   - Starting position in Cartesian coordinates.
%       X2   - Ending position in Cartesian coordinates.
%       RM   - Optional rotation matrix for coordinate transformation.
%
%   Output:
%       rho  - Radial distance in cylindrical coordinates.
%       theta - Azimuthal angle in cylindrical coordinates.
%       z    - Axial distance in cylindrical coordinates.
%

    % Set default rotation matrix to identity matrix
    R = eye(3);

    % Check if a rotation matrix is provided
    if nargin > 2
        R = RM;
    end

    % Calculate the difference in position vectors
    delta_r = X2 - X1;

    % Apply the rotation matrix to the position difference
    delta_r = R * delta_r';

    % Calculate cylindrical coordinates
    rho = sqrt(delta_r(1)^2 + delta_r(2)^2);
    theta = atan2(delta_r(2), delta_r(1));
    z = delta_r(3);

    % Convert theta to degrees if needed
    % theta_deg = rad2deg(theta);

    % Display the results (commented out in this version)
    % fprintf('Cylindrical Coordinates:\n');
    % fprintf('r = %.4f\n', rho);
    % fprintf('theta = %.4f degrees\n', theta_deg);
    % fprintf('z = %.4f\n', z);
end

