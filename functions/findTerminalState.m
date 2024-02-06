function [maxTimeLength, startTime, endTime, Thres] = findTerminalState(particle_part, RelError)
% FINDTERMINALSTATE identifies the terminal state of particles based on relative error.
%   [maxTimeLength, startTime, endTime, Thres] = FINDTERMINALSTATE(particle_part, RelError)
%   analyzes particle data to find the terminal state based on specified relative error.
%
%   Input:
%       particle_part - Structure array containing particle information.
%       RelError      - Relative error threshold for identifying the terminal state.
%
%   Output:
%       maxTimeLength - Maximum time length of the identified terminal state.
%       startTime     - Start time of the terminal state.
%       endTime       - End time of the terminal state.
%       Thres         - Relative error threshold used for identification.
%

%     % Maximum number of frames
%     Nframemax = 1e6;
% 
%     % Calculate experiment numbers for each particle
%     for i = 1:numel(particle_part)
%         nexp(i,:) = fix(particle_part(i).Tf/Nframemax) + 1;
%     end

%     % Unique experiment numbers
%     uni_exp = unique(nexp);

    % Call the findTimeRange function to get time range and threshold
    [maxTimeLength, startTime, endTime, Thres] = findTimeRange(particle_part.Ay, RelError);

%     % Adjust start and end times based on experiment numbers
%     startTime = startTime + Nframemax * (uni_exp - 1);
%     endTime = endTime + Nframemax * (uni_exp - 1);
end
