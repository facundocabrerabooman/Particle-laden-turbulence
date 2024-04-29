function splitVelCyl = slipVeloCylindrical_fcb(particle_part, tracers_part, neighborIdxAll, startTime, endTime, k)
% SLIPVELOCYLINDERICAL computes slip velocities in cylindrical coordinates
% for tracers near particles between specified start and end times.
%   slipVeloCylind = SLIPVELOCYLINDERICAL(particle_part, tracers_part,
%   neighborIdxAll, startTime, endTime, k) calculates slip velocities based
%   on input parameters and stores the results in the structure 'slipVeloCylind'.
%
%   Input:
%       particle_part   - Structure containing information about particles.
%       tracers_part    - Structure containing information about tracers.
%       neighborIdxAll  - Structure containing neighbor information.
%       startTime       - Start time for slip velocity calculation.
%       endTime         - End time for slip velocity calculation.
%       k               - Frame-related parameter for calculations.
%
%   Output:
%       slipVeloCylind  - Structure containing slip velocity information.
%           .rho        - Radial distance in cylindrical coordinates.
%           .theta      - Azimuthal angle in cylindrical coordinates.
%           .z          - Axial distance in cylindrical coordinates.
%           .Urel       - Relative slip velocity in cylindrical coordinates.
%

Nframemax = 1e7;

% Set default value for frame-related parameter
kexp = 1;
if nargin > 5
    kexp = k;
end

% Calculate start and end times within the frame limit
% tstart = mod(startTime(kexp), Nframemax);
% tend = mod(endTime(kexp), Nframemax);
tstart = 1;
tend = numel(neighborIdxAll);

for t = tstart:tend

    % Prepare indices for searching
    idx_front_back = sort(vertcat(neighborIdxAll(t).idx));

    % Initialize structure to store slip velocity information
    slipVeloCylind = struct('rho',[], 'theta',[], 'z',[], 'Urel',[]);

    % Iterate over indices/tracers
    for i = 1:numel(idx_front_back)
        idxp = particle_part.Tf == tracers_part.Tf(idx_front_back(i));
        idxt = idx_front_back(i);


        % Extract position and velocity information for particles and tracers
        Xp = [particle_part.Xf(idxp), particle_part.Yf(idxp), particle_part.Zf(idxp)];
        Vp = [particle_part.Vx(idxp), particle_part.Vy(idxp), particle_part.Vz(idxp)];
        Xf = [tracers_part.Xf(idxt), tracers_part.Yf(idxt), tracers_part.Zf(idxt)];
        Uf = [tracers_part.Vx(idxt), tracers_part.Vy(idxt), tracers_part.Vz(idxt)];

        % Rotate coordinates to align the particle velocity with the z-axis
        R = rotationMatrix(Vp, [0,0,1]);
        rXp = (R*Xp')';
        rXf = (R*Xf')';
        rVp = (R*Vp')';
        rUf = (R*Uf')';
        %         norm_rVp = rVp/norm(rVp);

        % Convert Cartesian coordinates to cylindrical coordinates
        [slipVeloCylind(i,:).rho, slipVeloCylind(i,:).theta, slipVeloCylind(i,:).z] = cartToCylind(rXp, rXf);

        % Calculate relative slip velocity in cylindrical coordinates
        slipVeloCylind(i,:).Urel = -(rUf - rVp); % FCB changed sign

        % Record time
        slipVeloCylind(i).t = t;

        splitVelCyl{t} = slipVeloCylind; 
    end
end
end
