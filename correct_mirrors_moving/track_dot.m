function [centers, radii] = track_dot(im, plotflag)

%%% Thresholds
maxAreaThreshold = 1500;
minAreaThreshold = 900;


imageHeight = size(im, 1);
imageWidth = size(im, 2);

mask = zeros(imageHeight, imageWidth, 'uint16');
mask(19:93, 1029:1117) = 1;
im_b = im .* mask;

% [centers, radii] = imfindcircles(imbinarize(im_b), [15, 40],'ObjectPolarity','dark')
%
% % Plot the circles on the original image
% imshow(imbinarize(im_b)); % Show the original image
% viscircles(centers, radii,'EdgeColor','b'); % Overlay the circles on the image

binaryImage = imcomplement(imbinarize(im_b));
cc = bwconncomp(binaryImage);
stats = regionprops(cc, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength');

filteredRegions = [];
for i = 1:length(stats)
    if stats(i).Area <= maxAreaThreshold && stats(i).Area >= minAreaThreshold
        filteredRegions = [filteredRegions; stats(i)];
    end
end

centroid = filteredRegions(1).Centroid;
radii = (filteredRegions(1).MajorAxisLength + filteredRegions(1).MinorAxisLength)/4;


if plotflag == 1
    figure(1);clf
    imshow(binaryImage);
    hold on;
    for i = 1:length(filteredRegions)
        centroid = filteredRegions(i).Centroid;
        majorAxisLength = filteredRegions(i).MajorAxisLength / 2;
        minorAxisLength = filteredRegions(i).MinorAxisLength / 2;
        plot(centroid(1), centroid(2), 'r*'); % Mark centroid
        rectangle('Position',[centroid(1)-majorAxisLength, centroid(2)-minorAxisLength, ...
            2*majorAxisLength, 2*minorAxisLength], 'EdgeColor', 'r', 'Curvature', [1 1]); % Draw ellipse
    end
    hold off;
end


end