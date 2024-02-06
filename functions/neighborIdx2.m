function [neighborIdxAll] = neighborIdx2(particle_part, tracer_part, idx1, idx2, Rmin, Rmax)
% NEIGHBORIDX2 computes neighboring tracers for a given particle in a specific frame.
%   [neighborIdxAll] = NEIGHBORIDX2(particle_part, tracer_part, idx1, idx2, Rmin, Rmax)
%   calculates neighboring tracers based on input parameters and stores the results
%   in the structure 'neighborIdxAll'.
%
%   Input:
%       particle_part - Structure containing information about particles.
%       tracer_part   - Structure containing information about tracers.
%       idx1          - Index of the current particle in 'particle_part'.
%       idx2          - Indices of tracers in the current frame in 'tracer_part'.
%       Rmin, Rmax    - Minimum and maximum distance constraints for neighbor search.
%
%   Output:
%       neighborIdxAll - Structure containing neighbor information.
%           .idx      - Indices of neighboring tracers within the specified distance range.
%           .d        - Distances between the current particle and neighboring tracers.
%           .Rmin     - Minimum distance constraint used in the search.
%           .Rmax     - Maximum distance constraint used in the search.
%           .idxfront - Indices of tracers in the front of the particle based on velocity.
%           .idxback  - Indices of tracers in the back of the particle based on velocity.
%

    % Extract information of the current particle
    particle.xf = particle_part.Xf(idx1);
    particle.yf = particle_part.Yf(idx1);
    particle.zf = particle_part.Zf(idx1);
    particle.vx = particle_part.Vx(idx1);
    particle.vy = particle_part.Vy(idx1);
    particle.vz = particle_part.Vz(idx1);

    % Extract information of tracers in the current frame
    tracer.xf = tracer_part.Xf(idx2);
    tracer.yf = tracer_part.Yf(idx2);
    tracer.zf = tracer_part.Zf(idx2);

    % Initialize structure to store neighbor information
    neighborIdxAll = struct('idx', [], 'd', []);
    neighborIdxAll.Rmin = Rmin;
    neighborIdxAll.Rmax = Rmax;

    % Calculate distances between the particle and each tracer
    for i = 1:numel(tracer.xf)
        d(i,:) = sqrt((particle.xf - tracer.xf(i))^2 + (particle.yf - tracer.yf(i))^2 + (particle.zf - tracer.zf(i))^2);
    end

    % Global neighboring tracers within the specified distance range
    idx01 = find(d > Rmin & d < Rmax);
    neighborIdxAll.idx = idx2(idx01);
    neighborIdxAll.d = d(idx01);

    % Front and back neighboring tracers based on velocity direction
    rpt.x = tracer.xf(idx01) - particle.xf;
    rpt.y = tracer.yf(idx01) - particle.yf;
    rpt.z = tracer.zf(idx01) - particle.zf;

    rv = particle.vx * rpt.x + particle.vy * rpt.y + particle.vz * rpt.z;
    neighborIdxAll.idxfront = neighborIdxAll.idx(rv > 0);
    neighborIdxAll.idxback = neighborIdxAll.idx(rv < 0);
end




