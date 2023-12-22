function mypolarpcolor(theta, radius, dataMatrix)
    % Custom function for polar pcolor-like plot

    % Create a polar plot with data points
    polar_ax = polaraxes;
    polarplot(theta(:), radius(:), 'k.'); % Plot polar points
    hold on;

    % Convert polar coordinates to Cartesian coordinates
    x = radius .* cos(theta);
    y = radius .* sin(theta);

    % Reshape data for pcolor
    numTheta = numel(unique(theta));
    numRadius = numel(unique(radius));
    reshapedData = reshape(dataMatrix, numRadius, numTheta);

    % Use pcolor on Cartesian coordinates
    ax = gca;
    cartesian_ax = axes('Position',ax.Position);
    pcolor(cartesian_ax, x, y, reshapedData);
    shading(cartesian_ax, 'interp'); % Interpolate colors for a smoother appearance
    colormap(cartesian_ax, jet);

    hold off;
end
