function [maxTimeLength, startTime, endTime,errorThreshold] = findTimeRange(signal, relativeErrorThreshold)
    % Input:
    % - signal: Time series signal
    % - relativeErrorThreshold: Maximum relative standard deviation allowed

    % Get the mean of the signal (using the mean of the end of the signal)
    endMeanSignal = mean(signal(1:end));  % Adjust the range as needed

    % Calculate the absolute error threshold based on the mean
    errorThreshold = abs(relativeErrorThreshold * endMeanSignal);

    % Initialize variables
    maxTimeLength = 0;
    startTime = 0;
    endTime = 0;

    % Get the length of the signal
    signalLength = length(signal);

    % Initialize variables for the current range
    currentStart = 1;
    currentEnd = 1;

    % Initialize standard deviation for the current range
    currentStdDev = std(signal(currentStart:currentEnd));

    while currentEnd <= signalLength
        % Check if the standard deviation is smaller than the threshold
        if currentStdDev < errorThreshold
            % Update if the current range is longer than the previous maximum
            if (currentEnd - currentStart + 1) > maxTimeLength
                maxTimeLength = currentEnd - currentStart + 1;
                startTime = currentStart;
                endTime = currentEnd;
            end

            % Expand the current range to the right
            currentEnd = currentEnd + 1;
            if currentEnd <= signalLength
                currentStdDev = std(signal(currentStart:currentEnd));
            end
        else
            % Shrink the current range from the left
            currentStart = currentStart + 1;
            if currentStart <= currentEnd
                currentStdDev = std(signal(currentStart:currentEnd));
            else
                % Expand the range if it has collapsed
                currentEnd = currentEnd + 1;
                if currentEnd <= signalLength
                    currentStdDev = std(signal(currentStart:currentEnd));
                end
            end
        end
    end

    % Display the results
%     disp(['Maximum Time Length: ', num2str(maxTimeLength)]);
%     disp(['Start Time: ', num2str(startTime)]);
%     disp(['End Time: ', num2str(endTime)]);
end
