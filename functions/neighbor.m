function neighborAll = neighbor(particle_part, tracer_part, Rmin, Rmax, kexp)
% NEIGHBOR computes neighbors for particles in a specific frame range.
%   neighborAll = NEIGHBOR(particle_part, tracer_part, Rmin, Rmax, kexp)
%   calculates neighbors based on input parameters and stores the results in
%   arrays 'neighborAll'.
%
%   Input:
%       particle_part - Information about particles.
%       tracer_part   - Information about tracers.
%       Rmin, Rmax    - Distance constraints for neighbor search.
%       kexp          - Frame-related parameter for calculations.
%
%   Output:
%       neighborAll   - Array containing all neighbor indices for each particle.

    % Maximum number of frames
    %Nframemax = 1e6;

    % Find indices of particles within the specified frame range
    %idxp = find(particle_part.Tf < kexp * Nframemax & particle_part.Tf > (kexp - 1) * Nframemax);
    idxp = 1:numel(particle_part.Tf);

    % Iterate over particles within the specified frame range
    for i = 1:numel(idxp)
        idx1 = idxp(i);

        % Find indices of tracers in the current frame
        idx2 = find(tracer_part.Tf == particle_part.Tf(idx1));

        % Calculate neighbors using a helper function
        neighborAll(i,:) = neighborIdx2(particle_part, tracer_part, idx1, idx2, Rmin, Rmax);
    end
end

