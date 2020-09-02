%% Create nice image plotting

load

% Firstly create closed loop, similar problem as bee FOV
% Just based on angles to start with, can add TSP for more complicated shapes
sortedEdgeSubscripts = pointsUnwrapped(edgeIndexList, :); 

% Mirror or claculations on full set of subscripts so they match sorted egde.
offSetFullSubscripts = pointsUnwrapped(:, :);

offSetFullSubscripts(:,1) = offSetFullSubscripts(:,1) - mean(sortedEdgeSubscripts(:,1));

offSetFullSubscripts(:,2) = offSetFullSubscripts(:,2) - mean(sortedEdgeSubscripts(:,2));

offSetSparseSubscripts = sparsePointsUnwrapped;

offSetSparseSubscripts(:,1) = offSetSparseSubscripts(:,1) - mean(sortedEdgeSubscripts(:,1));

offSetSparseSubscripts(:,2) = offSetSparseSubscripts(:,2) - mean(sortedEdgeSubscripts(:,2));

sortedEdgeSubscripts(:,1) = sortedEdgeSubscripts(:,1) - mean(sortedEdgeSubscripts(:,1));

sortedEdgeSubscripts(:,2) = sortedEdgeSubscripts(:,2) - mean(sortedEdgeSubscripts(:,2));

% Take sort edge list based on angle around origin 
    % Sequential angle sorting, fails on loop backs
    % Need to do before TSP other wise that can't cope
angles = atan2(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2));
        
[~,  tempIndex] = sort(angles);

% Shortest path solution, works really well!
[~, edgeIndexListSorted] = createSortedLoopwithTSP(double(sortedEdgeSubscripts(tempIndex,:)));

% Add first point at end to close.
edgeIndexListSorted = tempIndex([edgeIndexListSorted' edgeIndexListSorted(1)]);

sortedEdgeSubscripts = sortedEdgeSubscripts(edgeIndexListSorted, :);

% Transform crease flag as well.
tempSortedIndex = zeros(length(edgeIndexListSorted),1);

tempSortedIndex(creasePointsInEdge) = 1;

tempSortedIndex = find(tempSortedIndex(edgeIndexListSorted));

creaseIndexLeft = find(sortedEdgeSubscripts(tempSortedIndex,1) < 0);

creaseIndexLeft = tempSortedIndex(creaseIndexLeft);

creaseIndexRight = find(sortedEdgeSubscripts(tempSortedIndex,1) > 0);

creaseIndexRight = tempSortedIndex(creaseIndexRight);

figure; hold on;

plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2))

plot(sortedEdgeSubscripts(creaseIndexLeft,1), sortedEdgeSubscripts(creaseIndexLeft,2), 'bx-')

plot(sortedEdgeSubscripts(creaseIndexRight,1), sortedEdgeSubscripts(creaseIndexRight,2), 'rx-')

% Create 2D map for plotting.
xRange = ceil(max(sortedEdgeSubscripts(:,1)) - min(sortedEdgeSubscripts(:,1)) + 10);

% Mirror Y for image orientation 
sortedEdgeSubscripts(:,2) = -sortedEdgeSubscripts(:,2);

offSetFullSubscripts(:,2) = -offSetFullSubscripts(:,2);

offSetSparseSubscripts(:,2) = -offSetSparseSubscripts(:,2);

yRange = ceil(max(sortedEdgeSubscripts(:,2)) - min(sortedEdgeSubscripts(:,2)) + 10);

% Offset all matrices
offSetFullSubscripts(:,1) = offSetFullSubscripts(:,1) - min(sortedEdgeSubscripts(:,1)) + 5;

offSetFullSubscripts(:,2) = offSetFullSubscripts(:,2) - min(sortedEdgeSubscripts(:,2)) + 5;

offSetSparseSubscripts(:,1) = offSetSparseSubscripts(:,1) - min(sortedEdgeSubscripts(:,1)) + 5;

offSetSparseSubscripts(:,2) = offSetSparseSubscripts(:,2) - min(sortedEdgeSubscripts(:,2)) + 5;

sortedEdgeSubscripts(:,1) = sortedEdgeSubscripts(:,1) - min(sortedEdgeSubscripts(:,1)) + 5;

sortedEdgeSubscripts(:,2) = sortedEdgeSubscripts(:,2) - min(sortedEdgeSubscripts(:,2)) + 5;

% Create mask image.
image2Plot = zeros(xRange , yRange, 'logical');

% Test which points on map are within aleurone border.
[X, Y] = meshgrid(1:xRange, 1:yRange);

X = X(:); Y = Y(:);

[in, on] = inpolygon(X, Y, sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2));

inAleurone = in + on;

for iPixel = 1:length(in)
    
    if inAleurone(iPixel)
        
        image2Plot(X(iPixel), Y(iPixel)) = 1;
    end
end

% Do opening on image to tidy up edge points...
image2Plot = imopen(image2Plot, strel('disk', 2));

% Get index and subscripts again
inMap = find(image2Plot(:));

[XPointsIn, YPointsIn] = ind2sub([xRange, yRange], inMap);

%%% Should also identify which subscripts are properly in map (they can extend beyond border limits)
% Just use those in following

% Plot to test
figure; imshow(image2Plot'); hold on;

plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2), 'b', 'linewidth', 2)

plot(sortedEdgeSubscripts(creaseIndexLeft,1), sortedEdgeSubscripts(creaseIndexLeft,2), 'm*-')

plot(sortedEdgeSubscripts(creaseIndexRight,1), sortedEdgeSubscripts(creaseIndexRight,2), 'm*-');

plot(offSetFullSubscripts(:,1), offSetFullSubscripts(:,2), 'r.')

%plot(offSetSparseSubscripts(:,1), offSetSparseSubscripts(:,2), 'g.')

%% Get voxels ID for aleurone and endosperm

IDInterpolant = scatteredInterpolant( double([offSetFullSubscripts(:,1)' offSetSparseSubscripts(:,1)']'), ...
    double([offSetFullSubscripts(:,2)' offSetSparseSubscripts(:,2)']'),...
    [interpolatedIdentity' sparseIdentity']',...
    'nearest','nearest');

% Interpolate into map - keep unit8 to force integers
IDImage = zeros(xRange , yRange, 'uint8');

IDImage(inMap) = uint8(IDInterpolant(XPointsIn, YPointsIn)) + 1;

% Close a bit to smooth, but should avoid closing over holes
IDImage = imclose(IDImage, strel('disk',3));

%tempIDImage = imclose(IDImage, strel('disk',5)) - IDImage;
%figure; imshow(tempIDImage*100)

% Get specific points in each
inAleurone = find(IDImage(inMap) == 2);

XPointsInAleurone = XPointsIn(inAleurone); YPointsInAleurone = YPointsIn(inAleurone); 

inAleurone = inMap(inAleurone);

inEndosperm = find(IDImage(inMap) == 1);

XPointsInEndosperm = XPointsIn(inEndosperm); YPointsIEndosperm = YPointsIn(inEndosperm); 

inEndosperm = inMap(inEndosperm);

figure; imshow((IDImage')*100); hold on;

%´Take borders - endosperm to all
tempIm = IDImage;

% Add air then remove endosperm
tempIm(IDImage == 0) = 2;

tempIm(IDImage == 1) = 0;

tempIm = imdilate(tempIm, strel('disk',1));

endoBorderInds = find(IDImage == 1 & tempIm);

[endoBorderX, endoBorderY] = ind2sub([xRange , yRange], endoBorderInds);

hold on
plot(endoBorderX, endoBorderY, 'g.')

% Take borders - aleurone to air
tempIm = IDImage*0;

% Just add air
tempIm(IDImage == 0) = 1;

tempIm = imdilate(tempIm, strel('disk',1));

aleuroneBorderInds = find(IDImage == 2 & tempIm);

[aleuroneBorderX, aleuroneBorderY] = ind2sub([xRange , yRange], aleuroneBorderInds);

plot(aleuroneBorderX, aleuroneBorderY, 'y.')

plot(sortedEdgeSubscripts(creaseIndexLeft,1), sortedEdgeSubscripts(creaseIndexLeft,2), 'm-', 'linewidth', 2)

plot(sortedEdgeSubscripts(creaseIndexRight,1), sortedEdgeSubscripts(creaseIndexRight,2), 'm-', 'linewidth', 2);

%% Calculate thickness and intesntiy image
%figure; subplot(1,2,1); hist(thicknessByPoint,100)
%subplot(1,2,2); hist(averageIntensityByPoint,100)

valuesToUse = find(~isnan(thicknessByPoint));

valuesToSparse = find(~isnan(thicknessForSparse));

plotRefPoints = 0;
plotRefLines = 0;
refLines = {[1286, 2116], [1503, 1801]; ...
            [573.5 2077], [750.9 1745]; ...
            [322.2, 1009],[438.4 1654]; ...
            [2051 1106], [2128 635.3]};

figure; 
subplot(1,3,1);
[n,x] = hist([thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']',100);
plot(x,smooth(n)/sum(n),'linewidth',2)
xlabel('Thickness (um)')
ylabel('Percentage');

subplot(1,3,2);
[n,x] = hist([averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']',100);
plot(x,smooth(n)/sum(n),'linewidth',2)
xlabel('Intensity (AU)')

subplot(1,3,3);
plot([thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']', ...
    [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']', '.');

rValue = corr([thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']', ...
    [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']');

title(sprintf('r = %.2f', rValue))
xlabel('Thickness (um)'); ylabel('Intensity (AU)');


thicknessInterpolant = scatteredInterpolant( double([offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']'), ...
    double([offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']'),...
    [thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']',...
    'linear','nearest');

thicknessImage = zeros(xRange , yRange);

thicknessImage(inAleurone) = thicknessInterpolant(XPointsInAleurone, YPointsInAleurone);

warning('Colour range setting is not automated')

[max(thicknessByPoint) max(thicknessForSparse)]

thicknessImage = (thicknessImage-0)/(100-0);

figure; 
imshow(thicknessImage')
colormap(gray(100))
title('Thickness (um)'); hcb = colorbar; set(hcb,'Ticks', [0 0.5 1], 'TickLabels', {'0','40','80'})

hold on
plot(endoBorderX, endoBorderY, 'g.')
plot(aleuroneBorderX, aleuroneBorderY, 'y.')

if plotRefPoints
    plot(offSetFullSubscripts(:,1), offSetFullSubscripts(:,2), 'r.')
end

if plotRefLines
   nLines = size(refLines,1);
   
   for iLine = 1:nLines
       line1 = refLines{iLine,1};
       
       line2 = refLines{iLine,2};
       
       plot([line1(1) line2(1)], [line1(2) line2(2)], 'm')
   end
end

% Do for intensity.
intensityInterpolant = scatteredInterpolant( double([offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']'), ...
    double([offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']'),...
    [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']',...
    'linear','nearest');

intensityImage = zeros(xRange , yRange);

intensityImage(inAleurone) = intensityInterpolant(XPointsInAleurone, YPointsInAleurone);


[max(averageIntensityByPoint) max(averageIntensityForSparse)]
[min(averageIntensityByPoint) min(averageIntensityForSparse)]

intensityImage = (intensityImage-100)/(250-100);

figure;
imshow(intensityImage')
colormap(gray(100))
title('Intensity (AU)'); hcb = colorbar; set(hcb,'Ticks', [0 0.5 1], 'TickLabels', {'100','160','220'})

hold on
plot(endoBorderX, endoBorderY, 'g.')
plot(aleuroneBorderX, aleuroneBorderY, 'y.')

if plotRefPoints
    plot(offSetFullSubscripts(:,1), offSetFullSubscripts(:,2), 'r.')
end

if plotRefLines
   nLines = size(refLines,1);
   
   for iLine = 1:nLines
       line1 = refLines{iLine,1};
       
       line2 = refLines{iLine,2};
       
       plot([line1(1) line2(1)], [line1(2) line2(2)], 'm')
   end
end
%% Calculate intensity depth interpolant
% First pick dominant ID for each point - also set coords for interp
pointProfileIDMax = zeros(nPoints, numberOfBlocks)*NaN; 

pointXDepthCoord = zeros(nPoints, numberOfBlocks)*NaN; 
pointYDepthCoord = zeros(nPoints, numberOfBlocks)*NaN; 
pointZDepthCoord = zeros(nPoints, numberOfBlocks)*NaN; 

sparseProfileIDMax = zeros(nSparsePoints, numberOfBlocks)*NaN;

sparseXDepthCoord = zeros(nSparsePoints, numberOfBlocks)*NaN; 
sparseYDepthCoord = zeros(nSparsePoints, numberOfBlocks)*NaN; 
sparseZDepthCoord = zeros(nSparsePoints, numberOfBlocks)*NaN; 

for iBlock = 1:numberOfBlocks
    
    for jPoint = 1:nPoints
       % Simply take index of max as ID 
       [~, pointProfileIDMax(jPoint,iBlock) ] = max(pointProfileID(jPoint,iBlock,:));
       
       pointXDepthCoord(jPoint,iBlock) = offSetFullSubscripts(jPoint,1);
       pointYDepthCoord(jPoint,iBlock) = offSetFullSubscripts(jPoint,2);
       pointZDepthCoord(jPoint,iBlock) = (iBlock-1)*blockThickness;
    end
    
    for jPoint = 1:nSparsePoints
       % Simply take index of max as ID 
       [~, sparseProfileIDMax(jPoint,iBlock) ] = max(sparseProfileID(jPoint,iBlock,:));
       
       sparseXDepthCoord(jPoint,iBlock) = offSetSparseSubscripts(jPoint,1);
       sparseYDepthCoord(jPoint,iBlock) = offSetSparseSubscripts(jPoint,2);
       sparseZDepthCoord(jPoint,iBlock) = (iBlock-1)*blockThickness;
    end
end

% Subtract one as index offset 1 up from usual
pointProfileIDMax = pointProfileIDMax - 1;

sparseProfileIDMax = sparseProfileIDMax - 1;

% Take usable inds;
valuesToUse = find(~isnan(pointProfileIDMax));

valuesToSparse = find(~isnan(sparseProfileIDMax));

% Unwrap all
pointXDepthCoord = pointXDepthCoord(valuesToUse); sparseXDepthCoord = sparseXDepthCoord(valuesToSparse);

pointYDepthCoord = pointYDepthCoord(valuesToUse); sparseYDepthCoord = sparseYDepthCoord(valuesToSparse);

pointZDepthCoord = pointZDepthCoord(valuesToUse); sparseZDepthCoord = sparseZDepthCoord(valuesToSparse);

pointProfileIDMax = pointProfileIDMax(valuesToUse); sparseProfileIDMax = sparseProfileIDMax(valuesToSparse);

% Also unwrap intensity values for later
pointIntensityProfileUW = pointIntensityProfile(valuesToUse); sparseIntensityProfileUW = sparseIntensityProfile(valuesToSparse);

% Set up interpolant
depthIDInterpolant = scatteredInterpolant( double([pointXDepthCoord' sparseXDepthCoord']'), ...
    double([pointYDepthCoord' sparseYDepthCoord']'),...
    double([pointZDepthCoord' sparseZDepthCoord']'),...
    [pointProfileIDMax' sparseProfileIDMax']', 'nearest','nearest');

figure; hist([pointProfileIDMax' sparseProfileIDMax']')

% Interpolate by layer
IDProfileInterp = zeros(length(XPointsIn), numberOfBlocks)*NaN;

for iBlock = 1:numberOfBlocks
    IDProfileInterp(:, iBlock) = depthIDInterpolant(XPointsIn, YPointsIn, ...
        ones(length(XPointsIn),1)*(iBlock-1)*blockThickness);
end

%%
figure; 

colsID = [[0 0 0]; [1 0 1]; [1 1 0]; [0 1 0]; [0 0 1]];

% Resort cols ID
colsID = colsID([1:2, [ALEURONE_INDEX, ENDOSPERM_INDEX, GERM_INDEX]+2],:);

c = 1;
for iBlock = 1:numberOfBlocks/10:numberOfBlocks   
    subplot(numPlotRows,5,c); c = c+1;
    
    tempMap = zeros(xRange , yRange, 'uint8');
    
    tempMap(inMap) = IDProfileInterp(:, iBlock)+1;
    
    tempMap(1:4) = 1:4;
    
    colormap(colsID)
    
    imagesc(tempMap')
    axis off; axis equal
end
%% interpolate intensity profile in 3D and display depth maps

[max(pointIntensityProfile(:)) max(sparseIntensityProfile(:))]
[min(pointIntensityProfile(:)) min(sparseIntensityProfile(:))]

% Firstly make interpolant from non-nan values
valuesToUse = find(~isnan(pointIntensityProfileUW));

valuesToSparse = find(~isnan(sparseIntensityProfileUW));

intensityProfileInterpolant = scatteredInterpolant( double([pointXDepthCoord(valuesToUse)' sparseXDepthCoord(valuesToSparse)']'), ...
    double([pointYDepthCoord(valuesToUse)' sparseYDepthCoord(valuesToSparse)']'),...
    double([pointZDepthCoord(valuesToUse)' sparseZDepthCoord(valuesToSparse)']'),...
    [pointIntensityProfileUW(valuesToUse)' sparseIntensityProfileUW(valuesToSparse)']', 'linear','nearest');

% Interpolate by layer
intensityProfileInterp = zeros(length(XPointsIn), numberOfBlocks)*NaN;

for iBlock = 1:numberOfBlocks
    indsToGet = find(IDProfileInterp(:, iBlock) == ENDOSPERM_INDEX);
    
    intensityProfileInterp(indsToGet, iBlock) = intensityProfileInterpolant(XPointsIn(indsToGet), ...
        YPointsIn(indsToGet), ones(length(indsToGet),1)*(iBlock-1)*blockThickness);
end

%%
figure;

% Put into arrays and plot
c = 1;
for iBlock = 1:numberOfBlocks/10:numberOfBlocks 
   
    tempImage = zeros(xRange , yRange);
    
    % Have to solve w/ color image because of scaling
    tempImageCol = zeros(xRange , yRange, 3);
    
    indsToUse = find(IDProfileInterp(:, iBlock) == ENDOSPERM_INDEX);
    
    % Get points and wrap to range
    tempImage(inMap(indsToUse)) = (intensityProfileInterp(indsToUse, iBlock)-50)/(150-50);
    
    tempImageR = tempImage; tempImageG = tempImage; tempImageB = tempImage;
    
    % Get other IDs to add
    indsToUse = find(IDProfileInterp(:, iBlock) ~= ENDOSPERM_INDEX);
    
    % Note, offset of 2 required for this colormap indexing
    tempImageR(inMap(indsToUse)) = colsID(IDProfileInterp(indsToUse, iBlock)+2, 1);
    tempImageG(inMap(indsToUse)) = colsID(IDProfileInterp(indsToUse, iBlock)+2, 2);
    tempImageB(inMap(indsToUse)) = colsID(IDProfileInterp(indsToUse, iBlock)+2, 3);
    
    tempImageCol(:,:,1) = tempImageR;
    tempImageCol(:,:,2) = tempImageG;
    tempImageCol(:,:,3) = tempImageB;
    
    subplot(numPlotRows,5,c); c = c + 1;

    imshow( permute(tempImageCol, [2 1 3])); hold on;
    
    plot(endoBorderX, endoBorderY, 'g.')
    plot(aleuroneBorderX, aleuroneBorderY, 'y.')
    
    title(sprintf('%i - %i um', [(iBlock-1) iBlock]*blockThickness*VOXEL_SIZE));
end
    
figure
colormap(gray(100))
hcb = colorbar; set(hcb,'Ticks', [0 0.5 1], 'TickLabels', {'50','100','150'})
title('Profile intensity (AU)')
%% Look at intensity correlation between aleurone and depth layers

% Also plot histogram
figure;
hold on

cols = winterCoolColours;

for iBlock = 1:numberOfBlocks
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock))); 
    
    [n, x] = hist([pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)']', 100);
    
    plot(x, n/sum(n), 'color', cols(iBlock,:))
end
ylabel('Percentage');
xlabel('Intensity (AU)')

figure;
c = 1;
for iBlock = 1:numberOfBlocks/10:numberOfBlocks 
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)) & ~isnan(thicknessByPoint));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock)) & ~isnan(thicknessForSparse));
    
    alIntensity = [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']';
    
    blockIntensity = [pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)']';
    
    subplot(numPlotRows,5,c); c = c + 1;
    plot(alIntensity, blockIntensity, '.');
    
    rValue = corr(alIntensity, blockIntensity);
    title(sprintf('%i - %i, r = %.2f', [(iBlock-1) iBlock]*blockThickness*VOXEL_SIZE, rValue));
end

subplot(numPlotRows,5,6)
xlabel('Aleurone intensity (AU)')
ylabel('Endosperm intensity (AU)')

figure;
c = 1;
for iBlock = 1:numberOfBlocks1:numberOfBlocks/10:numberOfBlocks 
    alThickness = [thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']';
    
    blockIntensity = [pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)']';
    
    subplot(numPlotRows,5,c); c = c+1;
    plot(alThickness, blockIntensity, '.');
    
    rValue = corr(alThickness, blockIntensity);
    title(sprintf('%i - %i, r = %.2f', [(iBlock-1) iBlock]*blockThickness*VOXEL_SIZE, rValue));
end

subplot(numPlotRows,5,6)
xlabel('Aleurone thickness (um)')
ylabel('Endosperm intensity (AU)')

%% Test plot a slice
% Go between top, bottom reference points.
[lineX, lineY] = bresenham(offSetFullSubscripts(indexAbove,1), offSetFullSubscripts(indexAbove,2), ...
    offSetFullSubscripts(indexBelow,1), offSetFullSubscripts(indexBelow,2));

lineX2 = lineX; lineY2 = lineY;

% Try centre of crease indexes
% [lineX, lineY] = bresenham(sortedEdgeSubscripts(round(mean(creaseIndexLeft)),1), ...
%     sortedEdgeSubscripts(round(mean(creaseIndexLeft)),2), ...
%     sortedEdgeSubscripts(round(mean(creaseIndexRight)),1), ...
%     sortedEdgeSubscripts(round(mean(creaseIndexRight)),2));

% On test slice - get line between extemes (should be crease)
    % Should really sort and chain together rather than just taking between ends...
[~, minInd] = min(offSetSparseSubscripts(sparsePointsInSlice,1));
[~, maxInd] = max(offSetSparseSubscripts(sparsePointsInSlice,1));
[lineX, lineY] = bresenham(offSetSparseSubscripts(sparsePointsInSlice(minInd),1), offSetSparseSubscripts(sparsePointsInSlice(minInd),2), ...
    offSetSparseSubscripts(sparsePointsInSlice(maxInd),1), offSetSparseSubscripts(sparsePointsInSlice(maxInd),2));

lineDistance = sqrt((lineY - lineY(1)).^2 + (lineX-lineX(1)).^2);

% Get indexes into interpolation map.
[pointsToUse, indsOrig] = intersect([lineX, lineY], [XPointsIn, YPointsIn], 'rows');

% Checke line distance are sorted
[lineDistance, tempInd] = sort(lineDistance(indsOrig));

% Apply sort to others
pointsToUse = double(pointsToUse(tempInd,:));

% Check line distance step is consistent
if any(diff(lineDistance) > 1.3)
    figure; plot(diff(lineDistance))
    
   error('Line stepping over un-even distances') 
elseif any(diff(lineDistance) > 1)
   warning('Line stepping over un-even distances') 
end

% Just to test line
warning('Plot this on better map for reference')
figure; imshow(IDImage'*100); hold on;

% plot(offSetFullSubscripts(indexAbove,1), offSetFullSubscripts(indexAbove,2), 'cd')
% 
% plot(offSetFullSubscripts(indexBelow,1), offSetFullSubscripts(indexBelow,2), 'cd')

% plot(offSetSparseSubscripts(sparsePointsInSlice,1), ...
%     offSetSparseSubscripts(sparsePointsInSlice,2), '.')

plot(lineX, lineY, 'm')

plot(lineX2, lineY2, 'r')

% Load values into slice, and also aleurone thickness
profileSlice = zeros(length(lineDistance), numberOfBlocks)*NaN;

IDSlice = zeros(length(lineDistance), numberOfBlocks)*NaN;

for iBlock = 1:numberOfBlocks
    IDSlice(:, iBlock) = depthIDInterpolant(pointsToUse(:,1), pointsToUse(:,2),...
        ones(length(lineDistance),1)*(iBlock-1)*blockThickness); 
    
    indsToUse = find(IDSlice(:, iBlock) == ENDOSPERM_INDEX);
    
    profileSlice(indsToUse, iBlock) = intensityProfileInterpolant(pointsToUse(indsToUse,1), pointsToUse(indsToUse,2),...
        ones(length(indsToUse),1)*(iBlock-1)*blockThickness);  
end

thicknessSlice = thicknessInterpolant(pointsToUse(:,1), pointsToUse(:,2)); 

intensitySlice = intensityInterpolant(pointsToUse(:,1), pointsToUse(:,2));

profileSlice = (profileSlice-40)/(120-40);

% make colour image with other regions
tempSliceCol = zeros(length(lineDistance), numberOfBlocks*blockThickness, 3);

tempSliceR = profileSlice;  tempSliceG = profileSlice;  tempSliceB = profileSlice; 

indsToUse = find(IDSlice ~= ENDOSPERM_INDEX);

tempSliceR(indsToUse) = colsID(IDSlice(indsToUse)+2, 1);
tempSliceG(indsToUse) = colsID(IDSlice(indsToUse)+2, 2);
tempSliceB(indsToUse) = colsID(IDSlice(indsToUse)+2, 3);

tempSliceCol(:,:,1) = imresize(tempSliceR, [length(lineDistance), numberOfBlocks*blockThickness]);
tempSliceCol(:,:,2) = imresize(tempSliceG, [length(lineDistance), numberOfBlocks*blockThickness]);
tempSliceCol(:,:,3) = imresize(tempSliceB, [length(lineDistance), numberOfBlocks*blockThickness]);

% plotting
figure;

imshow( permute(tempSliceCol, [2 1 3]))

hold on; axis equal; set(gca, 'Clipping', 'off'); axis on

% Note, scale thickness to voxel size as remainder in voxels
area(lineDistance, -thicknessSlice/VOXEL_SIZE, 0, 'FaceColor', 'y');
xlim([0 length(lineDistance)])
ylim([-25 numberOfBlocks*blockThickness])

yTickVals = [-25, 0, 25, 50, 75];
set(gca, 'YTick', yTickVals, 'YTickLabel', -yTickVals*VOXEL_SIZE)
xTickVals = (1:10)/VOXEL_SIZE*1000;
set(gca, 'XTick', xTickVals, 'XTickLabel', xTickVals*VOXEL_SIZE/1000)

title('Slice into grain')
xlabel('Distance along slice (mm)')
ylabel('Distance down slice (um)')

%% Plot slices reference from map
%%% Points are specified in aleurone map section
nLines = size(refLines,1);
  
figTest = figure; hold on; axis equal
plot3(subscriptsToInterpolate(:,1), subscriptsToInterpolate(:,2), subscriptsToInterpolate(:,3), '.')

for iLine = 1:nLines
    figure(figTest)
    
   line1 = refLines{iLine,1};

   line2 = refLines{iLine,2};

   % Match each point
   [~, ind1] = min(sqrt((line1(1) - offSetFullSubscripts(:,1)).^2 + (line1(2) - offSetFullSubscripts(:,2)).^2));
   
   [~, ind2] = min(sqrt((line2(1) - offSetFullSubscripts(:,1)).^2 + (line2(2) - offSetFullSubscripts(:,2)).^2));
   
   % Plot each point and normal in 3D to test
   %%% May need to add a catch in case normal is empty later on
   normal1 = normalByPoint(ind1,:);
    
   point1 = subscriptsToInterpolate(ind1,:);

   plot3(point1(1), point1(2), point1(3), 'x')    
       
  line([0 50]*normal1(1) + point1(1),[0 50]*normal1(2) + point1(2), ...
        [0 50]*normal1(3) + point1(3)); 
    
   normal2 = normalByPoint(ind2,:);
    
   point2 = subscriptsToInterpolate(ind2,:);

   plot3(point2(1), point2(2), point2(3), 'x')    
       
   line([0 50]*normal2(1) + point2(1), [0 50]*normal2(2) + point2(2), ...
        [0 50]*normal2(3) + point2(3)); 
    
   plot3([point1(1) point2(1)], [point1(2) point2(2)], [point1(3) point2(3)], 'linewidth',2)
   
   pointLine = point1 - point2;
   
   pointDist = norm(pointLine);
   
   pointLine = pointLine/norm(pointLine);

   normalLIne = (normal1 + normal2)/2;
   
   % Take rotation matrix to align x-axis
   rMat = matrix2rotatevectors(pointLine, [1 0 0]);
   % Rotate y-axis
   planeLine = [0 1 0]*rMat;
   
   %Solve for best rotation angle
   angles = -pi:pi/360:pi;
   
   angleError = zeros(length(angles), 1);
   
   for jAngle = 1:length(angles)
      % Rotate b around a by theta,then take angle to c
      testLine = planeLine*vrrotvec2mat([pointLine angles(jAngle)]);
      
      angleError(jAngle) = atan2(norm(cross(testLine,normalLIne)),...
          dot(testLine,normalLIne));
   end
   
   % Apply rotation
   [~, minI] = min(angleError);
   
   planeLine = planeLine*vrrotvec2mat([pointLine angles(minI)]);
   
   line([0 100]*planeLine(1) + point2(1), [0 100]*planeLine(2) + point2(2), ...
        [0 100]*planeLine(3) + point2(3), 'color', 'm'); 
   
   %%% Should centre plane points
   
   testX = -200:(round(pointDist)+200);
   testY = -100:100;
   
   planeX = zeros(length(testX), length(testY));
   planeY = zeros(length(testX), length(testY));
   planeZ = zeros(length(testX), length(testY));
   
   for jPoint = 1:length(testX)
      planeX(jPoint,:) = testX(jPoint)*pointLine(1) + testY(:)*planeLine(1) + point2(1);
      
      planeY(jPoint,:) = testX(jPoint)*pointLine(2) + testY(:)*planeLine(2) + point2(2);
      
      planeZ(jPoint,:) = testX(jPoint)*pointLine(3) + testY(:)*planeLine(3) + point2(3);
   end
   
   
   plot3(planeX(:), planeY(:), planeZ(:), '.')
   
   interpIm = interp3(single(greyVolumeAligned),planeY,planeX,planeZ,'linear');
   figure;
   imshow(uint8(interpIm'))
end
   
%% Profile plotting
%%% May wish to add option of normalizing to radius
nBins = 20;

% Try for aleurone thickenss and intensity to start with.

% Points should be roughly equally distributed
areaPerPoint = aleuroneExteriorArea/(sum(interpolatedIdentity) + sum(sparseIdentity))/10^6;
% 10mm long 2.5 mm diamter tube area is about 75 mm2, seems good ball park

% For verticalthickness and intensity
valuesToUse = find(~isnan(thicknessByPoint));

valuesToSparse = find(~isnan(thicknessForSparse));

[verticalThicknessMean, verticalCoords, ~, verticalCount, vertBins] = calculateProfile([thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']',...
    -[offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']', nBins, []);

[verticalIntensityMean] = calculateProfile([averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']',...
    -[offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']', vertBins, []);

% Expected area error
(length(valuesToUse)+length(valuesToSparse))/(sum(interpolatedIdentity) + sum(sparseIdentity))
% Actual area error
sum(verticalCount)/(sum(interpolatedIdentity) + sum(sparseIdentity))

vertFigure = figure; 

vertOffset = 0; -min(verticalCoords)+max(verticalCoords);

subplot(1,3,1)
plot(verticalThicknessMean, (verticalCoords+vertOffset)*VOXEL_SIZE/1000);
xlabel('Average aleurone thickness (um)')
ylabel('Distance from base (mm)')
xlim([0 60]); ylim([0 12])

subplot(1,3,2)
plot(verticalIntensityMean, (verticalCoords+vertOffset)*VOXEL_SIZE/1000);
xlabel('Average aleurone intensity (AU)')
title('Vertical profiles')
xlim([0 200]); ylim([0 12])

subplot(1,3,3)
plot(verticalCount*areaPerPoint, (verticalCoords+vertOffset)*VOXEL_SIZE/1000);
xlabel('Average aleurone area (mm2)')
xlim([0 10]); ylim([0 12])

% For horizonalt thickness and intensity

xCent = mean(offSetFullSubscripts([indexAbove indexBelow],1));

% Area is large at bottom becasue of points clustering on base
[horizontalThicknessMean, horizontalCoords, ~, horizontalCount, hoirzBins] = calculateProfile([thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']',...
    [offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']', nBins, xCent);

%[subscriptsToInterpolate(valuesToUse,3)' subscriptsForSparse(valuesToSparse,3)']

[horizontalIntensityMean] = calculateProfile([averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']',...
    [offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']', hoirzBins, xCent);

figure;
subplot(1,3,1);
plot(horizontalCoords*VOXEL_SIZE/1000, horizontalThicknessMean);
ylabel('Average aleurone thickness (um)')
ylim([0 60]); xlim([-5 5])

subplot(1,3,2)
plot(horizontalCoords*VOXEL_SIZE/1000, horizontalIntensityMean);
ylabel('Average aleurone intensity (AU)')
title('Horizontal profiles')
xlabel('Distance from centreline (mm)')
ylim([0 200]); xlim([-5 5])

subplot(1,3,3)
plot(horizontalCoords*VOXEL_SIZE/1000, horizontalCount*areaPerPoint);
ylabel('Average aleurone area (mm2)')
xlim([-5 5]); ylim([0 5])

% plot vertical nad horizontal profiles for endosperm intensity

cols = winterCoolColours;

%%% Could also plot intensity profile in aleurone for comparison to exposed endosperm

figure;
for iBlock = 1:numberOfBlocks
   % Mask to area under aleurone 
   valuesToUseEndo = find(~isnan(pointIntensityProfile(:,iBlock)) & ~isnan(thicknessByPoint));

   valuesToSparseEndo = find(~isnan(sparseIntensityProfile(:,iBlock)) & ~isnan(thicknessForSparse));
   
    % Calc certical, using bins from al
    [verticalIntensityMean] = calculateProfile([pointIntensityProfile(valuesToUseEndo, iBlock)' sparseIntensityProfile(valuesToSparseEndo, iBlock)']',...
        -[offSetFullSubscripts(valuesToUseEndo,2)' offSetSparseSubscripts(valuesToSparseEndo,2)']', vertBins, []);
    
    subplot(2,2,1); hold on
    plot(verticalIntensityMean, (verticalCoords+vertOffset)*VOXEL_SIZE/1000, 'color', cols(iBlock,:));

    [horizontalIntensityMean] = calculateProfile([pointIntensityProfile(valuesToUseEndo, iBlock)' sparseIntensityProfile(valuesToSparseEndo, iBlock)']',...
        [offSetFullSubscripts(valuesToUseEndo,1)' offSetSparseSubscripts(valuesToSparseEndo,1)']', hoirzBins, xCent);

    subplot(2,2,2); hold on
    plot(horizontalCoords*VOXEL_SIZE/1000, horizontalIntensityMean, 'color', cols(iBlock,:));
    
    % Take full area
    valuesToUseEndo = find(~isnan(pointIntensityProfile(:,iBlock)) & isnan(thicknessByPoint));

    valuesToSparseEndo = find(~isnan(sparseIntensityProfile(:,iBlock)) & isnan(thicknessForSparse));
   
    % Calc certical, using bins from al
    [verticalIntensityMean] = calculateProfile([pointIntensityProfile(valuesToUseEndo, iBlock)' sparseIntensityProfile(valuesToSparseEndo, iBlock)']',...
        -[offSetFullSubscripts(valuesToUseEndo,2)' offSetSparseSubscripts(valuesToSparseEndo,2)']', vertBins, []);
    
    subplot(2,2,3); hold on
    plot(verticalIntensityMean, (verticalCoords+vertOffset)*VOXEL_SIZE/1000, 'color', cols(iBlock,:));

    [horizontalIntensityMean] = calculateProfile([pointIntensityProfile(valuesToUseEndo, iBlock)' sparseIntensityProfile(valuesToSparseEndo, iBlock)']',...
        [offSetFullSubscripts(valuesToUseEndo,1)' offSetSparseSubscripts(valuesToSparseEndo,1)']', hoirzBins, xCent);

    subplot(2,2,4); hold on
    plot(horizontalCoords*VOXEL_SIZE/1000, horizontalIntensityMean, 'color', cols(iBlock,:));
end

subplot(2,2,1); hold on
xlabel('Average endosperm intensity (AU)')
ylabel('Distance from base (mm)')
xlim([0 100]); ylim([0 12])
title('Vertical layers under aleurone')

subplot(2,2,2); hold on
ylabel('Average endosperm intensity (AU)')
xlabel('Distance from centreline (mm)')
ylim([0 100]); xlim([-5 5])
title('Horizontal layers under aleurone')

subplot(2,2,3); hold on
xlabel('Average endosperm intensity (AU)')
ylabel('Distance from base (mm)')
xlim([0 200]); ylim([0 12])
title('Vertical layers under exposed endosperm')

subplot(2,2,4); hold on
ylabel('Average endosperm intensity (AU)')
xlabel('Distance from centreline (mm)')
ylim([0 200]); xlim([-5 5])
title('Horizontal layers under exposed endosperm')
%% Calculate distance error images
% Need to distance calculate between points in 2D.
distanceMatrix2D = zeros(nPoints, nPoints)*NaN;

% Get 2D indexes of points
pointIndex2D = sub2ind([xRange, yRange], round(offSetFullSubscripts(:,1)), ...
    round(offSetFullSubscripts(:,2)));

% Not guaranteed that borders will be inside, so need to replace with index of closest point.
for iPoint = 1:nPoints
    
   if image2Plot(pointIndex2D(iPoint)) == 0 
        
       [~, minInd] = min( sqrt((offSetFullSubscripts(iPoint,1)-XPointsIn).^2 + ...
           (offSetFullSubscripts(iPoint,2) - YPointsIn).^2));
       
       pointIndex2D(iPoint) = sub2ind([xRange, yRange], XPointsIn(minInd), YPointsIn(minInd));
   end
end

% Get 2D distance maps.
tic
parfor iPoint = 1:nPoints
    dMap = bwdistgeodesic(image2Plot, pointIndex2D(iPoint), 'quasi-euclidean');
    
    distanceMatrix2D(iPoint, :) = dMap(pointIndex2D);
end
toc

figure; 
temp = abs(distanceMatrix(:) - distanceMatrix2D(:));
subplot(2,2,1); hist(temp(~isinf(temp)),100); title('all')

% Get min error, max error and difference

minErrorList = zeros(nPoints,1)*NaN;

meanErrorList = zeros(nPoints,1)*NaN;

maxErrorList = zeros(nPoints,1)*NaN;

for iPoint = 1:nPoints
   tempError =  abs(distanceMatrix(iPoint, :) - distanceMatrix2D(iPoint, :));
   
   tempError(iPoint) = [];
   
   if sum(~isinf(tempError))
       minErrorList(iPoint) = min(tempError(~isinf(tempError)));
       
       meanErrorList(iPoint) = mean(tempError(~isinf(tempError)));

       maxErrorList(iPoint) = max(tempError(~isinf(tempError)));
   end
end

subplot(2,2,3); hist(minErrorList,100); title('min')

subplot(2,2,4); hist(maxErrorList,100); title('max')

subplot(2,2,2); hist(meanErrorList,100); title('Mean')



% Create error images for mean.
errorInterpolant = scatteredInterpolant( double(offSetFullSubscripts(:,1)), ...
   double(offSetFullSubscripts(:,2)), meanErrorList, 'nearest','nearest');

errorImage = zeros(xRange , yRange);

errorImage(inMap) = errorInterpolant(XPointsIn, YPointsIn);

figure; subplot(1,3,1); imshow(errorImage'/round(max(meanErrorList))); title('mean')

% Create error images for min.
errorInterpolant = scatteredInterpolant( double(offSetFullSubscripts(:,1)), ...
    double(offSetFullSubscripts(:,2)), minErrorList, 'nearest','nearest');

errorImage = zeros(xRange , yRange);

errorImage(inMap) = errorInterpolant(XPointsIn, YPointsIn);

subplot(1,3,2); imshow(errorImage'/round(max(minErrorList))); title('min')

% Create error images for max.
errorInterpolant = scatteredInterpolant( double(offSetFullSubscripts(:,1)), ...
    double(offSetFullSubscripts(:,2)), maxErrorList, 'nearest','nearest');

errorImage = zeros(xRange , yRange);

errorImage(inMap) = errorInterpolant(XPointsIn, YPointsIn);

subplot(1,3,3); imshow(errorImage'/round(max(maxErrorList))); title('max')

%%% Could also plot where maximum location of error goes to, if far away
%%% it's probably not important.
