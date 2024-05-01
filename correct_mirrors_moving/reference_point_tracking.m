clc

figure(1);clf
% Define the dimensions of your image
imageHeight = size(x9_c1, 1);
imageWidth = size(x9_c1, 2);

% Create a mask with ones inside the specified rectangle and zeros outside
mask = zeros(imageHeight, imageWidth, 'uint16'); % Use the same data type as your image
mask(19:93, 1029:1117) = 1;

% Apply the mask to your image
x9_c1_b = x9_c1 .* mask;

% [centers, radii] = imfindcircles(imbinarize(x9_c1_b), [15, 40],'ObjectPolarity','dark')
% 
% % Plot the circles on the original image
% imshow(imbinarize(x9_c1_b)); % Show the original image
% viscircles(centers, radii,'EdgeColor','b'); % Overlay the circles on the image

binaryImage = imcomplement(imbinarize(x9_c1_b));
cc = bwconncomp(binaryImage);
stats = regionprops(cc, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength');

maxAreaThreshold = 1500; 
minAreaThreshold = 900; 
filteredRegions = [];
for i = 1:length(stats)
    if stats(i).Area <= maxAreaThreshold && stats(i).Area >= minAreaThreshold
        filteredRegions = [filteredRegions; stats(i)];
    end
end

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


%% 
% Define the dimensions of your image
imageHeight = size(x13_c1, 1);
imageWidth = size(x13_c1, 2);

% Create a mask with ones inside the specified rectangle and zeros outside
mask = zeros(imageHeight, imageWidth, 'uint16'); % Use the same data type as your image
mask(19:93, 1029:1117) = 1;

% Apply the mask to your image
x13_c1_b = x13_c1 .* mask;

[centers, radii] = imfindcircles(imbinarize(x13_c1_b), [15, 40],'ObjectPolarity','dark')

figure(2)
% Plot the circles on the original image
imshow(imbinarize(x13_c1_b)); % Show the original image
viscircles(centers, radii,'EdgeColor','r'); % Overlay the circles on the image


