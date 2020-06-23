% To do
    %Normalize intensity between multiple grains based on germ histogram

% Test load and plot of grain stack.
clc; clear; close all

% Other: Om_1_6
labelDirectory = 'C:\Grain_temp\Om1_6\Labels';
greyDirectory = 'C:\Grain_temp\Om1_6\Images';
dataName = 'OM16';
ALEURONE_INDEX = 2; % Material index.
ENDOSPERM_INDEX = 1;
GERM_INDEX = 3;
turnXY = -1; %Set to get correct rotation of crease, should face down
roundOnLoop = 16;
removeLowest = 1;
padLow = [20 0 0];
padHigh = [0 0 0];

%OB_6 test
%labelDirectory = 'C:\Grain_temp\OB6_test\Labels';
%greyDirectory = 'C:\Grain_temp\OB6_test\Images';
%dataName = 'OB6';
% ALEURONE_INDEX = 1; % Material index.
% ENDOSPERM_INDEX = 2;
% GERM_INDEX = 3;
%turnXY = 1;
%roundOnLoop = 64; % Can be set back to 8?
%removeLowest = 2;
%padLow = [0 0 0];
%padHigh = [0 0 0];

% Loading tiff stack takes a while.
grainVolume = loadtiffstack(labelDirectory, 1);

% pad volume
if any(padLow); grainVolume = padarray(grainVolume, padLow, 0, 'pre'); end

if any(padHigh); grainVolume = padarray(grainVolume, padHigh, 0, 'post'); end

volumeSize = size(grainVolume);

% Create stucturing elements for sepecific connections
STREL_6_CONNECTED = strel('sphere', 1); 
temp = STREL_6_CONNECTED.Neighborhood; 
temp(:,:,2) = 1; temp([1 3],2,:) = 1; temp(2,[1 3],:) = 1;
STREL_18_CONNECTED = strel('arbitrary', temp); 
STREL_26_CONNECTED = strel('cube', 3); 

% Called imopen in parts as error in Mex file.
s = settings; s.images.UseHalide.TemporaryValue = false;
grainVolume = imopen(grainVolume, STREL_6_CONNECTED);
s = settings; s.images.UseHalide.TemporaryValue = true;

VOXEL_SIZE = 4;

% Points on maps, in voxels
edgeDistance = 20; %10

surfaceDistance = 50; %50

%Not, sparse distance must be less then others for code to work
%OM16, will thrash on 2 during distance map, maybe remove parallel
sparsePointsDistance = 4; %3 

% Flags
testPlotNormals = 0;

%OM1_6 has a short cut, OB6 does not, probably best to stop... will add 10 minutes
stopShortCuts = 1;

maxAleuroneThickness = 75; %In voxels - previous limit at 20

blockThickness = 8/VOXEL_SIZE;

depthToCalculate = 120/VOXEL_SIZE;

% Try to make these all integer values given Voxel thickenss
numberOfBlocks = depthToCalculate/blockThickness;

numPlotRows = ceil(numberOfBlocks/5);

winterCoolColours = [(0:numberOfBlocks-1)', fliplr(0:numberOfBlocks-1)' (0:numberOfBlocks-1)']/(numberOfBlocks-1);

% Curve is calculated at this interval
curveStep = 5;

% Density of nrb spline interpolation
    % 10 - ~min, 5 - ~15 min
% Note that original interpolation is every xth slice, 
% so interp is at curveStepxstep slices
%%% Set to 5 for good runs
nrbStep = 10;

% offsetCut in by this amount, otherwise just steps in 10 from ends.
    %i.e. total cut in 10 + offsetCut
offsetCut = 40;
%% Get voxels on exterior of orignal grain by overlap to exterior.
% Get pixels in grain.
grainIndexList = find(grainVolume);

nIndex = length(grainIndexList(1:50:end)); grainSubscriptArray = zeros(nIndex, 3);

[grainSubscriptArray(:,1), grainSubscriptArray(:,2), grainSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, grainIndexList(1:50:end));

% Get long axis of grain with PCA.
grainAxisArray = pca(grainSubscriptArray);

grainLongAxis = grainAxisArray(:,1)';

grainCreaseAxis = grainAxisArray(:,3)';

grainCenter = mean(grainSubscriptArray);

transform2Vertical = matrix2rotatevectors([0, 0, 1], grainLongAxis);

transform2Up = matrix2rotatevectors([0, turnXY, 0], grainCreaseAxis*transform2Vertical);

% Apply transformation to grain in volume.
temp = zeros(4,4); temp(4,4) = 1;

M1 = make_transformation_matrix(grainCenter - volumeSize/2);

M2 = temp; M2(1:3,1:3) = transform2Vertical;

M3 = temp; M3(1:3,1:3) = transform2Up;

M4 = make_transformation_matrix(volumeSize/2 - grainCenter);

% Tried to test if bounds get exceeded - results do not match affine-transform...
% grainSubscriptArrayRotating = [grainSubscriptArray-volumeSize/2, ones(nIndex,1)];
% 
% grainSubscriptArrayRotating = grainSubscriptArrayRotating*M1*M2*M3*M4;
% 
% grainSubscriptArrayRotating = grainSubscriptArrayRotating(:,1:3) + volumeSize/2;

% Pad adjust centres
volumeSize = size(grainVolume);

grainVolumeAligned = uint8( affine_transform_full(single(grainVolume), M1*M2*M3*M4, 5));

%%% Note small holes appear after affine transform, close fixes these. 
%%% Probably bug with nearest neighbour interp I added to affine transform c file
grainVolumeAligned = imclose(grainVolumeAligned, STREL_6_CONNECTED);

% Check for touch on borders
if sum(sum( any(grainVolumeAligned(1,:,:)))); error('Touch on x 1'); end

if sum(sum( any(grainVolumeAligned(end,:,:)))); error('Touch on X end'); end

if sum(sum( any(grainVolumeAligned(:,1,:)))); error('Touch on Y 1'); end

if sum(sum( any(grainVolumeAligned(:,end,:)))); error('Touch on Y end'); end

if sum(sum( any(grainVolumeAligned(:,:,1)))); error('Touch on Z 1'); end

if sum(sum( any(grainVolumeAligned(:,:,end)))); error('Touch on Z end'); end

%  figure; imshow(grainVolumeAligned(:,:,1375)*100)

% Get exterior voxles
%%% ADD FUNCTION for this and replace in loop.
grainExterior = logical(~grainVolumeAligned);

% Exterior is outer volume, and grown into other volume.
grainExterior = imdilate(grainExterior, STREL_18_CONNECTED);

grainExterior = grainExterior & grainVolumeAligned;

% Take largest connected region.
tempCC = bwconncomp(grainExterior, 18);

tempStats = regionprops(tempCC, 'PixelIdxList');

% Get number of voxels in each region. 
nRegions = length(tempStats); voxelsPerRegionArray = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
end

% Largest will generally be much larger than others.
[~, tempIndex] = max(voxelsPerRegionArray); tempStats(tempIndex) = [];

% Remove other regions from volume.
for iRegion = 1:nRegions-1
    grainExterior(tempStats(iRegion).PixelIdxList) = 0;
end

% Get surface indices, then convert to subscripts.
surfaceIndexList = find(grainExterior);

nIndex = length(surfaceIndexList); grainSurfaceSubscriptArray = zeros(nIndex, 3);

[grainSurfaceSubscriptArray(:,1), grainSurfaceSubscriptArray(:,2), grainSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, surfaceIndexList);

% Get aleurone surface subscripts.
aleuroneSurfaceIndexList = find(grainExterior & grainVolumeAligned == ALEURONE_INDEX);

nIndex = length(aleuroneSurfaceIndexList); aleuroneSurfaceSubscriptArray = zeros(nIndex, 3);

[aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), aleuroneSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, aleuroneSurfaceIndexList);

clear grainIndexList, clear grainVolume

figure; 
% hold on; axis equal; set(gca, 'Clipping', 'off')
% plot3(grainSurfaceSubscriptArray(:,1), grainSurfaceSubscriptArray(:,2), grainSurfaceSubscriptArray(:,3), 'b.')
subplot(1,2,1);
imshow(grainVolumeAligned(:,:,500)*100); title('500')
subplot(1,2,2);
imshow(grainVolumeAligned(:,:,2000)*100); title('2000')
%% Fill beneath crease on each z-slices to get shaped.
% Get ranges to check.
zRange = [ceil( min(grainSurfaceSubscriptArray(:,3))) floor( max(grainSurfaceSubscriptArray(:,3)))];

zRangeAleurone = [ceil( min(aleuroneSurfaceSubscriptArray(:,3))) floor( max(aleuroneSurfaceSubscriptArray(:,3)))];

xRange = [ceil( min(grainSurfaceSubscriptArray(:,1))) floor( max(grainSurfaceSubscriptArray(:,1)))];

yRange = [ceil( min(grainSurfaceSubscriptArray(:,2))) floor( max(grainSurfaceSubscriptArray(:,2)))];

% Get aleurone distribution by slices.
aleuroneVoxelsBySlice = zeros(volumeSize(3),1);
for iSlice = zRangeAleurone(1):zRangeAleurone(2)
    
    aleuroneVoxelsBySlice(iSlice) = sum( sum(grainVolumeAligned(:,:,iSlice) == ALEURONE_INDEX));
    
end

% Set zRange for testing.
zRangeAleurone = find(aleuroneVoxelsBySlice); %[500 2000]; %[30 2110]

zRangeAleurone = [zRangeAleurone(1)+50 zRangeAleurone(end)-50];

figure; plot(aleuroneVoxelsBySlice); hold on
plot(zRangeAleurone, aleuroneVoxelsBySlice(zRangeAleurone), 'rx')

% Sum of zeros 1st, then max of zeros.
creaseProfileBySlice = zeros(volumeSize(3),volumeSize(1),2);

boundsLineBySlice = zeros(volumeSize(3),2)*NaN;

exteriorVolume = ones(volumeSize, 'logical');

% Identfiy split bounds by slice.
for iSlice = zRangeAleurone(1):zRangeAleurone(2)
    
    % Step along x in columns of y.
    for jColumn = 1:volumeSize(1)
        
        temp = find(grainVolumeAligned(jColumn,:,iSlice) == ALEURONE_INDEX | ...
                grainVolumeAligned(jColumn,:,iSlice) == ENDOSPERM_INDEX);
        
        if ~isempty(temp)
            % Find highest point.
            grainTop = max(temp);
            
            % Find zeros underneath.
            zerosIndexList = find(grainVolumeAligned(jColumn,1:grainTop,iSlice) == 0);

            exteriorVolume(jColumn,zerosIndexList,iSlice) = 0;

            % Save profile.
            creaseProfileBySlice(iSlice, jColumn, 1) = length(zerosIndexList);
            
            creaseProfileBySlice(iSlice, jColumn, 2) = max(zerosIndexList);
        end
    end

    % Find bulges on grain by getting peaks on slice. Invert and rectify.
    tempColumn = -1*creaseProfileBySlice(iSlice, :, 1);
    
    tempColumn(tempColumn < 0) = tempColumn(tempColumn < 0) - min(tempColumn(tempColumn < 0));

    % Suggest minimum width, works better with general at start.
    if iSlice < mean(zRangeAleurone)
        widthTemp = (xRange(2)-xRange(1))/3;
    else
        widthTemp = find(tempColumn);
        
        widthTemp = (widthTemp(end)-widthTemp(1))/3;
    end
    
    % Suggest max height.
    heightTemp = max([max(tempColumn)*0.25 50]);
    
    % Width and height minimums are somewhat arbitary, seems robust for OB6.
    [~, tempPeakIndexList] = findpeaks(tempColumn, 'MinPeakDistance', widthTemp, 'MinPeakHeight', heightTemp); 
    
    % Ignore if more than 2
    if length(tempPeakIndexList) == 2
        % Places bounds at average between center and edges.
        boundsLineBySlice(iSlice,:) = tempPeakIndexList;
    end
end

clear exteriorVolume

% Interp any missing values.
tempMissing = find( isnan(boundsLineBySlice(:,1)));

tempHave = find( ~isnan(boundsLineBySlice(:,1)));

boundsLineBySlice(tempMissing,1) = interp1(tempHave, boundsLineBySlice(tempHave,1), ...
    tempMissing, 'linear');

boundsLineBySlice(tempMissing,2) = interp1(tempHave, boundsLineBySlice(tempHave,2), ...
    tempMissing, 'linear');

% Smooth.
tempHave = find( ~isnan(boundsLineBySlice(:,1)));

boundsLineBySlice(tempHave,1) = round( smooth(boundsLineBySlice(tempHave,1)));

boundsLineBySlice(tempHave,2) = round( smooth(boundsLineBySlice(tempHave,2)));

figure; imshow(creaseProfileBySlice(:,:,1)/300); hold on

plot(boundsLineBySlice(:,1), 1:volumeSize(3), 'r.')

plot(boundsLineBySlice(:,2), 1:volumeSize(3), 'r.')

% Estimate centreline.
centreLine = zeros(volumeSize(3),1)*NaN;

for iSlice = zRangeAleurone(1):zRangeAleurone(2)
    
    if all( ~isnan( boundsLineBySlice(iSlice,:)))
        
        temp = boundsLineBySlice(iSlice,1):boundsLineBySlice(iSlice,2);
    
        centreLine(iSlice) = wmean(temp, creaseProfileBySlice(iSlice, temp, 1).^2);
        
    end
end

tempHave = find( ~isnan(centreLine));

centreLine(tempHave) = round(smooth(centreLine(tempHave)));

plot(centreLine, 1:volumeSize(3), 'g.')

% Get height along centreline.
centreLineHeightProfile = zeros(volumeSize(3),1)*NaN;

for iSlice = zRangeAleurone(1):zRangeAleurone(2)
    
    if all( ~isnan( boundsLineBySlice(iSlice,:)))
        
        centreLineHeightProfile(iSlice) = creaseProfileBySlice(iSlice,centreLine(iSlice),2);
        
    end
end
%% Create new volume to calculate loop under. 
% First set new bounds.
zBoundsNew = [zRange(1)-10 zRange(2)+10];

if zBoundsNew(1) < 1; zBoundsNew(1) = 1; end

if zBoundsNew(2) > volumeSize(3); zBoundsNew(2) = volumeSize(3); end

yBoundsNew = [yRange(1)-10 yRange(2)+10];
    
if yBoundsNew(1) < 1; yBoundsNew(1) = 1; end

if yBoundsNew(2) > volumeSize(2); yBoundsNew(2) = volumeSize(2); end

xBoundsNew = [min(boundsLineBySlice(:,1))-10 max(boundsLineBySlice(:,2))+10];

if xBoundsNew(1) < 1; xBoundsNew(1) = 1; end

if xBoundsNew(2) > volumeSize(1); xBoundsNew(2) = volumeSize(1); end

smallGrainVolume = grainVolumeAligned(xBoundsNew(1):xBoundsNew(2), ...
    yBoundsNew(1):yBoundsNew(2), zBoundsNew(1):zBoundsNew(2));

smallGrainExterior = grainExterior(xBoundsNew(1):xBoundsNew(2), ...
    yBoundsNew(1):yBoundsNew(2), zBoundsNew(1):zBoundsNew(2));

smallVolumeSize = size(smallGrainVolume);

%warning('Orignal volumes cleared')

%clear grainVolumeAligned, clear grainExterior



% Find loop start on top of grain. Firstly step in from top of grain
zTopOfLoop = zRange(1)-zBoundsNew(1)+1 + 10;

% Then take centre X position.
sliceIndexList = find(smallGrainVolume(:,:,zTopOfLoop));

nIndex = length(sliceIndexList); sliceSubscriptArray = zeros(nIndex, 3);

[sliceSubscriptArray(:,1), sliceSubscriptArray(:,2), sliceSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, sliceIndexList);

xTopOfLoop = round( mean(sliceSubscriptArray(:,1)));

% Then take lowest Y in line.
tempIndex = find(sliceSubscriptArray(:,1) == xTopOfLoop);

yTopOfLoop = min(sliceSubscriptArray(tempIndex,2))-1;

% Repeat for bottom of loop.
zBottomOfLoop = zRange(2)-zBoundsNew(1)+1 - 10;

sliceIndexList = find(smallGrainVolume(:,:,zBottomOfLoop));

nIndex = length(sliceIndexList); sliceSubscriptArray = zeros(nIndex, 3);

[sliceSubscriptArray(:,1), sliceSubscriptArray(:,2), sliceSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, sliceIndexList);

xBottomOfLoop = round( mean(sliceSubscriptArray(:,1)));

tempIndex = find(sliceSubscriptArray(:,1) == xBottomOfLoop);

yBottomOfLoop = min(sliceSubscriptArray(tempIndex,2))-1;
%% Calc distance map from plane underneath grain.
% Make mask image for grain.
grainMask = logical(smallGrainVolume);

[testX, testY] = meshgrid(1:smallVolumeSize(1), 1:smallVolumeSize(2));

% Step through and cut above centre line;
for iSlice = 1:smallVolumeSize(3)
    
    if ~isnan( centreLineHeightProfile(iSlice+zBoundsNew(1)-1))
        
        cutAbove = centreLineHeightProfile(iSlice+zBoundsNew(1)-1)-yBoundsNew(1)+1+50;
        
        if cutAbove < smallVolumeSize(2)
            
            grainMask(:, cutAbove:end, iSlice) = 1;
        end 
    else
        % Make polygon to block out triangle above centre of mass of pixels
        sliceIndexList = find(smallGrainVolume(:,:,iSlice));

        [xSubscriptList, ySubscriptList] = ind2sub(smallVolumeSize, sliceIndexList); 
        
        if ~isempty(ySubscriptList)
            % Get vertices in polygons.
            cutFrom = [min(xSubscriptList) max(xSubscriptList)];
           
            cutFrom = round( (cutFrom-mean(cutFrom))*0.5 + mean(cutFrom));
            
            cutAbove = round( mean(ySubscriptList));
            
            cutMiddle = (smallVolumeSize(2)-cutAbove)/3 + cutAbove;
            
            polyX = [cutFrom(1) cutFrom(2) smallVolumeSize(1) smallVolumeSize(1) 1 1 cutFrom(1)]';
            
            polyY = [cutAbove cutAbove cutMiddle smallVolumeSize(2) smallVolumeSize(2) cutMiddle cutAbove]';

            inPolyIndex = find( inpolygon(testX,testY,polyX,polyY));  
                        
            % Set indexes in poly
            if iSlice > zRange(1)-zBoundsNew(1)+1 && iSlice < zRange(2)-zBoundsNew(1)+1
                %If inside bounds, just set and specific slice

                inVolIndex = sub2ind(smallVolumeSize, testX(inPolyIndex), ...
                testY(inPolyIndex), ones(length(inPolyIndex),1)*iSlice);
            
                grainMask(inVolIndex) = 1;
            
            elseif iSlice == zRange(1)-zBoundsNew(1)+1
                % Fill to top of volume.
                for jStep = 1:iSlice
                    
                    inVolIndex = sub2ind(smallVolumeSize, testX(inPolyIndex), ...
                    testY(inPolyIndex), ones(length(inPolyIndex),1)*jStep);

                    grainMask(inVolIndex) = 1;
                end
                
            elseif iSlice == zRange(2)-zBoundsNew(1)+1
                % Fill to bottom of volume.
                for jStep = iSlice:smallVolumeSize(3)
                    
                    inVolIndex = sub2ind(smallVolumeSize, testX(inPolyIndex), ...
                    testY(inPolyIndex), ones(length(inPolyIndex),1)*jStep);

                    grainMask(inVolIndex) = 1;
                end
            end
        end
    end
end

% Calculate distance maps 
tic
distanceFromGrain = bwdist(grainMask, 'quasi-euclidean');
toc

% Need to set interior to max to prevent traversal. 
distanceFromGrain(grainMask) = Inf;

% Graydist takes a long time. Squared to force closer to grain (works?).
tic
distanceFromTop = graydist(distanceFromGrain.^2, sub2ind(smallVolumeSize, ...
    xTopOfLoop, yTopOfLoop, zTopOfLoop), 'quasi-euclidean');
toc

tic
distanceFromBottom = graydist(distanceFromGrain.^2, sub2ind(smallVolumeSize, ...
    xBottomOfLoop, yBottomOfLoop, zBottomOfLoop),'quasi-euclidean');
toc

dMap = distanceFromTop + distanceFromBottom;

% Round to lower precision to prevent floating point errors.
%%% Note, changed from 64 to 8 as loop doesn't reach end,
%%% Causes some floaters, but no gaps. Not guaranteed to work for all...
dMapRound = round(dMap * roundOnLoop)/roundOnLoop;

dMapRound(isnan(dMapRound)) = Inf;

% Reginal minima should define connection.
loopVolume = imregionalmin(dMapRound);

% Take largest connected region of loop.
tempCC = bwconncomp(loopVolume, 26);

tempStats = regionprops(tempCC, 'PixelIdxList');

% Get number of voxels in each region. 
nRegions = length(tempStats); voxelsPerRegionArray = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
end

% Largest will generally be much larger than others.
[~, tempIndex] = max(voxelsPerRegionArray); tempStats(tempIndex) = [];

% Remove other regions from volume.
for iRegion = 1:nRegions-1
    loopVolume(tempStats(iRegion).PixelIdxList) = 0;
end

loopIndexList = find(loopVolume);

nIndex = length(loopIndexList); loopSubscriptArray = zeros(nIndex, 3);

[loopSubscriptArray(:,1), loopSubscriptArray(:,2), loopSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, loopIndexList);

figure; subplot(1,2,1); hold on; axis equal

plot3(loopSubscriptArray(:,1), loopSubscriptArray(:,2), loopSubscriptArray(:,3), 'b.')

plot3(xTopOfLoop, yTopOfLoop, zTopOfLoop, 'rx')

plot3(xBottomOfLoop, yBottomOfLoop, zBottomOfLoop, 'gx')

subplot(1,2,2)

hist(dMap(loopIndexList),100)

% Check top and bottom of loop are in loop volume
if ~loopVolume(xTopOfLoop, yTopOfLoop, zTopOfLoop) | ... 
        ~loopVolume(xBottomOfLoop, yBottomOfLoop, zBottomOfLoop)
    
    % may need to adjsut round on loop value
    error('Loop did not reach ends, need to insepct')
end

clear dMapRound, clear distanceFromTop, clear distanceFromBottom
%% Can have a short cut if part of loop concave which is not solved in following section

if stopShortCuts
   % Take two largest values from histogram 
   [nValues, xCenter] = hist(dMap(loopIndexList),100);

   xDiff = diff(xCenter);
   
   % For some reasons values not equally spaced...
   xDiff = mean(xDiff (xDiff ~= 0));
   
   [nValuesSorted, valuesSortedIndex] = sort(nValues, 'descend');
   
   firstInds = find(dMap(loopIndexList) >= xCenter(valuesSortedIndex(1))& ...
       dMap(loopIndexList) < xCenter(valuesSortedIndex(1)) + xDiff);
   
   secondInds = find(dMap(loopIndexList) >= xCenter(valuesSortedIndex(2)) & ...
       dMap(loopIndexList) < xCenter(valuesSortedIndex(2)) + xDiff);
   
   %nValuesSorted(1:2)
   
   %[length(firstInds) length(secondInds)]
   
   %Select lowest to remove
   %%% Strictly speaking, the one to remove should have higher distance
   %%% values, but this may not be guaranteed
   
   if mean(loopSubscriptArray(firstInds,2)) < mean(loopSubscriptArray(secondInds,2))
   
       toRemove = firstInds;
   else
       toRemove = secondInds;
   end
   
   subplot(1,2,1)
   
   plot3(loopSubscriptArray(toRemove,1), loopSubscriptArray(toRemove,2), ...
       loopSubscriptArray(toRemove,3), 'rx')
   
   loopVolume(loopIndexList(toRemove)) = [];
   
   loopIndexList(toRemove) = [];
   
   loopSubscriptArray(toRemove, 1) = [];
   
   loopSubscriptArray(toRemove, 2) = [];
   
   loopSubscriptArray(toRemove, 3) = [];
end
clear dMap
%% Take top of connected curve.
%Provided both ends of loop are connect this kind of finds one path and
%ignores floaters...

coordinatesOfLoop = zeros(smallVolumeSize(3)*2,3)*NaN;

coordinatesOfLoop(1,:) = [xTopOfLoop, yTopOfLoop, zTopOfLoop];

[testX, testY, testZ] = meshgrid(-1:1, -1:1, -1:1);

testX(14) = []; testY(14) = []; testZ(14) = [];

roughXCentre = round((xTopOfLoop+xBottomOfLoop)/2);

tic
distanceFromEnd = bwdistgeodesic(loopVolume, sub2ind(smallVolumeSize, ...
  xBottomOfLoop, yBottomOfLoop, zBottomOfLoop), 'quasi-euclidean');
toc

% Step along loop towards end.
for iStep = 2:smallVolumeSize(3)*2
    % Get all distances to the adjacent points of previous point.
    distanceArray = zeros(26,1);
    
    distanceFromGrainArray = zeros(26,1);
    
    for jPoint = 1:26
        
        distanceArray(jPoint) = distanceFromEnd(coordinatesOfLoop(iStep-1,1)+testX(jPoint), ...
               coordinatesOfLoop(iStep-1,2)+testY(jPoint), coordinatesOfLoop(iStep-1,3)+testZ(jPoint)); 
           
        distanceFromGrainArray(jPoint) = distanceFromGrain(coordinatesOfLoop(iStep-1,1)+testX(jPoint), ...
               coordinatesOfLoop(iStep-1,2)+testY(jPoint), coordinatesOfLoop(iStep-1,3)+testZ(jPoint));   
    end

    % If just one point, take shortest.
    tempIndex = find(distanceArray == min(distanceArray));
    
    if length(tempIndex) == 1
        
        coordinatesOfLoop(iStep,:) = coordinatesOfLoop(iStep-1,:) + ...
            [testX(tempIndex) testY(tempIndex) testZ(tempIndex)];
        
        %plot3(coordinatesOfLoop(iStep,1), coordinatesOfLoop(iStep,2), coordinatesOfLoop(iStep,3), 'kx');
        
    elseif length(tempIndex) > 1
        % If more than 1 take, closest to grain
        %%% Add closest test as curve that diverged from main grain could
        %%% be slected, but didn't help much...
        
        tempIndex2 = find(distanceFromGrainArray(tempIndex) == min(distanceFromGrainArray(tempIndex)));
       
        if length(tempIndex2) == 1
            coordinatesOfLoop(iStep,:) = coordinatesOfLoop(iStep-1,:) + ...
                    [testX(tempIndex(tempIndex2)) testY(tempIndex(tempIndex2)) testZ(tempIndex(tempIndex2))];

        else
            % If still more than 1 point, take highest.
            tempIndex3 = find(testY(tempIndex(tempIndex2)) == max(testY(tempIndex(tempIndex2))));

            if length(tempIndex3) == 1

                coordinatesOfLoop(iStep,:) = coordinatesOfLoop(iStep-1,:) + ...
                    [testX(tempIndex(tempIndex2(tempIndex3))) testY(tempIndex(tempIndex2(tempIndex3))) ...
                    testZ(tempIndex(tempIndex2(tempIndex3)))];

                %plot3(coordinatesOfLoop(iStep,1), coordinatesOfLoop(iStep,2), coordinatesOfLoop(iStep,3), 'go');

            else

               % If still, still more than 1, take closest to centre. 
               tempDiffToCenter = abs(coordinatesOfLoop(iStep-1,1) + testX(tempIndex(tempIndex2(tempIndex3))) - roughXCentre);

               tempIndex4 = find(tempDiffToCenter == min(tempDiffToCenter));

               if length(tempIndex4) == 1
                   coordinatesOfLoop(iStep,:) = coordinatesOfLoop(iStep-1,:) + ...
                        [testX(tempIndex(tempIndex2(tempIndex3(tempIndex4)))) testY(tempIndex(tempIndex2(tempIndex3(tempIndex4)))) ...
                        testZ(tempIndex(tempIndex2(tempIndex3(tempIndex4))))];
               else
                   error('Multiple near centre points...');
               end
            end 
        end
    end
    
    % If end point reached, break early.
    if all((coordinatesOfLoop(iStep,:)-[xBottomOfLoop, yBottomOfLoop, zBottomOfLoop]) == 0)
        
        break;
        
    end
end

coordinatesOfLoop(isnan(coordinatesOfLoop(:,1)),:) = [];

clear distanceFromEnd

%Clear loop and add replace with thin values.
loopVolume(:) = 0;

tempIndexList = sub2ind(smallVolumeSize, coordinatesOfLoop(:,1), coordinatesOfLoop(:,2), coordinatesOfLoop(:,3));

loopVolume(tempIndexList) = 1;

loopIndexList = find(loopVolume);

nIndex = length(loopIndexList); loopSubscriptArray = zeros(nIndex, 3);

[loopSubscriptArray(:,1), loopSubscriptArray(:,2), loopSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, loopIndexList);

plot3(loopSubscriptArray(:,1), loopSubscriptArray(:,2), loopSubscriptArray(:,3), 'm.')

%% Now find centre-curve by slice.

% Transform distance from grain so centrelines are emphasized.
% Removed, seems to work better with linear distance
%distanceFromGrain = distanceFromGrain.^2;

distanceFromGrain(grainMask) = 0;

centreCurveVolume = zeros(smallVolumeSize, 'logical');

slicesToUse = [zTopOfLoop:curveStep:(zBottomOfLoop-curveStep+1) zBottomOfLoop];

if slicesToUse(end) ~= zBottomOfLoop
    slicesToUse = [slicesToUse zBottomOfLoop];
end

for iSlice = slicesToUse
   
    % Check if exisiting crease points of crease on slice.
    indexList = find(loopVolume(:,:,iSlice));
    
    if ~isempty(indexList)
        % Fast marching does a good job of finding good cut through grain.
        % But can jump a bit between slices if morphology shifts thinnest
        % crossing and sometimes is confused by diagonal openings.
    
        % Get points on loop.
        nIndex = length(indexList); tempSubscriptArray = zeros(nIndex, 2);

        [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(smallVolumeSize(1:2), indexList);

        % Just use highest point
        if nIndex > 1
            [~, maxInd] = max(tempSubscriptArray(:,2));

            tempSubscriptArray = tempSubscriptArray(maxInd,:);
        end

        % Calculate shortest path from centre of base plane to loop.
        distanceFromBase = msfm(double(distanceFromGrain(:,:,iSlice)*1000+0.001), [roughXCentre; 1]); 

        curveLine = shortestpath_robust(distanceFromBase,tempSubscriptArray',...
            [roughXCentre; 1], 0.5);

        %Use Bresenham algorithm to take all voxels touched by line.
        curveLineFull = zeros(size(curveLine,1)*4,2)*NaN;

        tempIndex = 1;

        for jStep = 1:length(curveLine)-1

            % Check distance is small.
            if sqrt((curveLine(jStep+1,1)-curveLine(jStep,1))^2 + ...
                    (curveLine(jStep+1,2)-curveLine(jStep,2))^2) < 5

                [tempX, tempY] = bresenham(curveLine(jStep,1),curveLine(jStep,2),...
                    curveLine(jStep+1,1),curveLine(jStep+1,2));

                curveLineFull(tempIndex:tempIndex+length(tempX)-1,:) = [tempX, tempY];

                tempIndex = tempIndex + length(tempX);
            else
               % Not expecting (really) large steps to occur
               error('Large distance!'); 
            end
        end

        curveLineFull(isnan(curveLineFull(:,1)), :) = [];

        curveLineFull = unique(curveLineFull, 'rows');

        %Voxelize line.
        curveSlice = zeros(smallVolumeSize(1:2));

        tempIndexList = sub2ind(smallVolumeSize(1:2), curveLineFull(:,1), curveLineFull(:,2));

        curveSlice(tempIndexList) = 1;

        tempMask = grainMask(:,:,iSlice);

        % Removed test for multiple connections to base plane, fast marching prevents.

        % Test that line links top and bottom.
        if curveSlice(roughXCentre,1) && ...
                curveSlice(tempSubscriptArray(1), tempSubscriptArray(2))
            % Removed thining and multi-part removal code, fast marching prevents.

            % Add slice of curve and loop to curve volume.
            centreCurveVolume(:,:,iSlice) = curveSlice | loopVolume(:,:,iSlice);

            %plot3(curveLineFull(:,1), curveLineFull(:,2), ...
            %    iSlice*ones(size(curveLineFull,1),1), 'k.')
            
            % Vertical continuity test removed, fast marching should fix.
        else
           % Start and end don't connect, just add exisiting loop. 
           % centreCurveVolume(:,:,iSlice) = loopVolume(:,:,iSlice);

           % Better not to ad anything.
        end
    end
end

clear distanceFromGrain, clear grainMask

%% Interpolate form using lofted b-spline.
curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, curveIndexList);

% Need to create array like meshgrid, but pad with nans because of variable length.
% Each row does not need to have same max Y value or same number of entries.

% Firstly get largest dimensions.
zValues = unique(curveSubscriptArray(:,3));

maxZLength = length( unique(curveSubscriptArray(:,3)));

maxYLength = 0;

for iSlice = 1:maxZLength
   indsInSlice = find(curveSubscriptArray(:,3) == zValues(iSlice));
   
   if maxYLength < length(indsInSlice)
       
       maxYLength = length(indsInSlice);
   end
end

% Then put curve subscripts into grid.
xArray = zeros(maxYLength, maxZLength)*NaN;

yArray = zeros(maxYLength, maxZLength)*NaN;

zArray = zeros(maxYLength, maxZLength)*NaN;

for iSlice = 1:maxZLength
    
    indsInSlice = find(curveSubscriptArray(:,3) == zValues(iSlice));
    
    xArray(1:length(indsInSlice), iSlice) = curveSubscriptArray(indsInSlice,1);
    
    yArray(1:length(indsInSlice), iSlice) = curveSubscriptArray(indsInSlice,2);
    
    zArray(1:length(indsInSlice), iSlice) = curveSubscriptArray(indsInSlice,3);
end

% Note - this is incredibly slow!
%%% This could be solved by doing sequentially for small sequences and then recombining.
%%% Do say, for 1-50 and 50-100 then combine. Should test 25-75 produces same result in overlap.
warning('Not using all points in nrbloft')

% Get subset to interpolat
toInterp = 1:nrbStep:size(xArray,2);

% Add final point if not included
if toInterp(end) ~= size(xArray,2)
    toInterp = [toInterp size(xArray,2)];
end

tic
curveNrb = nrbloft_nan(xArray(:, toInterp), yArray(:, toInterp), zArray(:, toInterp), 2);
toc

% 1st dimension is Y, 2nd dimension is Z
nrbplot(curveNrb, [10, 50]);
%% Put interpolated values into curve volume.
%zToInterp = zBottomOfLoop-zTopOfLoop+1;
zToInterp = max(max(zArray(:, toInterp)));

%yToInterp = max(loopSubscriptArray(:,2));
yToInterp = max(max(yArray(:, toInterp)));

% Multiply by two for effective half step, like fast marching.
p = nrbeval(curveNrb,{linspace(0.0,1.0,yToInterp*2) linspace(0.0,1.0,zToInterp*2)});
curveX = p(1,:,:); curveX = round(curveX(:));
curveY = p(2,:,:); curveY = round(curveY(:));
curveZ = p(3,:,:); curveZ = round(curveZ(:));

[~, indList] = unique([curveX, curveY, curveZ], 'rows');

curveX = curveX(indList); curveY = curveY(indList); curveZ = curveZ(indList);

% Find max Y for each Z value
maxYValues = zeros(smallVolumeSize(3),1);

for iSlice = zTopOfLoop:zToInterp %zBottomOfLoop
    
   inds = find(curveZ == iSlice);
   
   if ~isempty(inds)
       
       maxYValues(iSlice) = max(curveY(inds));
       
       % May be duplicate Y max points with different X values. 
       % Not a big problem as each will be filled up until the highest X in following loop.
       
   else
      error('No Y value') 
   end
end

% Fill up to X on each point. Is actually making a blockwe will get surface from
centreCurveVolume = zeros(smallVolumeSize, 'uint8');

for iPoint = 1:length(indList)
    
    centreCurveVolume(1:curveX(iPoint), curveY(iPoint), curveZ(iPoint)) = 1;
    
    %If max y for a given slice, fill up as well
    if curveY(iPoint) == maxYValues(curveZ(iPoint))
        
       for jRow = (maxYValues(curveZ(iPoint))+1):smallVolumeSize(2)
           
          centreCurveVolume(1:curveX(iPoint), jRow, curveZ(iPoint)) = 1;
           
       end
    end
end

% Copy ends along volume.
for iSlice = 1:zTopOfLoop
    
    centreCurveVolume(:,:,iSlice) = centreCurveVolume(:,:,zTopOfLoop);
end

for iSlice = zToInterp:smallVolumeSize(3)
    
    centreCurveVolume(:,:,iSlice) = centreCurveVolume(:,:,zToInterp);
end

% Take difference between volume erode and volume - should be continous surface
centreCurveVolume = uint8(centreCurveVolume - imerode(centreCurveVolume, STREL_26_CONNECTED));

% Remove ends from volume.
centreCurveVolume(:,:,1:(zTopOfLoop-1)) = 0;

centreCurveVolume(:,:,(zToInterp+1):smallVolumeSize(3)) = 0;

% Remove points above curve from volume
for iSlice = zTopOfLoop:zToInterp
    
   for jRow = (maxYValues(iSlice)+1):smallVolumeSize(2)

      centreCurveVolume(:, jRow, iSlice) = 0;

   end
end

% Get indexes again
curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, curveIndexList);

plot3(curveSubscriptArray(:,1), curveSubscriptArray(:,2), ...
    curveSubscriptArray(:,3), 'm.');
%% Extend curve, then place into main volume and re-calcualte exterior
%Extend curve to include front and end of grain.

for iSlice = 1:zTopOfLoop
    %Copy end of loop if border volume is not empty
    if ~isempty( find( smallGrainExterior(:,:,iSlice), 1))
        
        centreCurveVolume(:,:,iSlice) = centreCurveVolume(:,:,zTopOfLoop);
    end
end

for iSlice = zToInterp:smallVolumeSize(3)
    if ~isempty( find( smallGrainExterior(:,:,iSlice), 1))
        
        centreCurveVolume(:,:,iSlice) = centreCurveVolume(:,:,zToInterp);
    end
end

% Get subscripts for curve 
curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, curveIndexList);

% Place curve cut back into full volume
% Correct subscripts back to original coordinates
curveSubscriptArray(:,1) = curveSubscriptArray(:,1) + xBoundsNew(1) - 1;

curveSubscriptArray(:,2) = curveSubscriptArray(:,2) + yBoundsNew(1) - 1;

curveSubscriptArray(:,3) = curveSubscriptArray(:,3) + zBoundsNew(1) - 1;

curveIndexList = sub2ind(volumeSize, curveSubscriptArray(:,1), curveSubscriptArray(:,2), ...
    curveSubscriptArray(:,3));

refInds = find(grainVolumeAligned(curveIndexList));

[~, ~, tempSlice] = ind2sub(volumeSize, refInds);

% Cut label volume along base curve
grainVolumeAligned(curveIndexList) = 0;

% Recalculate exterior - if not, won't be closed along cut crease al 
grainExterior = logical(~grainVolumeAligned);

grainExterior = imdilate(grainExterior, STREL_18_CONNECTED);

grainExterior = grainExterior & grainVolumeAligned;

tempCC = bwconncomp(grainExterior, 18);

tempStats = regionprops(tempCC, 'PixelIdxList');

% Get number of voxels in each region. 
nRegions = length(tempStats); voxelsPerRegionArray = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
end

% Largest will generally be much larger than others.
[~, tempIndex] = max(voxelsPerRegionArray); tempStats(tempIndex) = [];

% Remove other regions from volume.
for iRegion = 1:nRegions-1
    grainExterior(tempStats(iRegion).PixelIdxList) = 0;
end

%% Now fill up form curve into exterior 
% First, resize curve volume
centreCurveVolume = zeros(volumeSize, 'uint8');

centreCurveVolume(curveIndexList) = 1;

% Idea was for all of curve to be under grain as full loop is, and then this cuts up
% But, can be that curve moved in x, intersects border, and doesnt get cut up
%%% Could test down from top as well to re-ID exisiting  grain intersects
for iSlice = 1:volumeSize(3)
    % Find y positions of curve on slice.   
    tempIndexList = find(centreCurveVolume(:,:,iSlice));

    if ~isempty(tempIndexList)
        nIndex = length(tempIndexList); tempSubscriptArray = zeros(nIndex, 2);

        [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(volumeSize(1:2), tempIndexList);
        
        [~, maxIndexList] = max(tempSubscriptArray(:,2));
        
        maxIndexList = find(tempSubscriptArray(:,2) == tempSubscriptArray(maxIndexList,2));
        
        % If multiple max y values on one slice cut up from each to ensure connectivity.
        
        for jIndex = 1:length(maxIndexList)
            % Get values in grain and edge volumes.
            volumeColumn = grainVolumeAligned(tempSubscriptArray(maxIndexList(jIndex),1), ...
                (tempSubscriptArray(maxIndexList(jIndex),2)+1):end, iSlice);

            exteriorColumn = grainExterior(tempSubscriptArray(maxIndexList(jIndex),1), ...
                (tempSubscriptArray(maxIndexList(jIndex),2)+1):end, iSlice);

            %Step up column and add into curve volume if air or edge
            for kPoint = 1:length(volumeColumn)

                if volumeColumn(kPoint) == 0 || exteriorColumn(kPoint) 
                    %|| volumeColumn(kPoint) == ALEURONE_INDEX
                
                    % Set value depending on region    
                    if iSlice > zTopOfLoop && iSlice < zToInterp    
                        centreCurveVolume(tempSubscriptArray(maxIndexList(jIndex),1), ...
                            tempSubscriptArray(maxIndexList(jIndex),2)+kPoint, iSlice) = 2;
                    else
                        centreCurveVolume(tempSubscriptArray(maxIndexList(jIndex),1), ...
                        tempSubscriptArray(maxIndexList(jIndex),2)+kPoint, iSlice) = 3;
                    end
%                     plot3(tempSubscriptArray(maxIndexList(jIndex),1),...
%                         tempSubscriptArray(maxIndexList(jIndex),2)+kPoint, iSlice, 'k.')

                % Stop once non- edge or air reach to prevent disconnection.
                % Only do if in main portion of grain, borders should extend up
                elseif iSlice > (zTopOfLoop + offsetCut) && iSlice < (zToInterp - offsetCut)
                    %if volumeColumn(kPoint) == ENDOSPERM_INDEX || volumeColumn(kPoint) == GERM_INDEX
                    
                    % plot3(tempSubscriptArray(maxIndexList(jIndex),1), tempSubscriptArray(maxIndexList(jIndex),2)+kPoint-1, iSlice, 'mx')

                   break
                end
            end
        end
    end
end

% Get indexes again
curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, curveIndexList);

mainCurveInds = find(centreCurveVolume(curveIndexList) == 1);

topCurveInds = find(centreCurveVolume(curveIndexList) == 2);

endCurveInds = find(centreCurveVolume(curveIndexList) == 3);

mainInVolInds = find(grainVolumeAligned(curveIndexList));

%% Test if curve ends are connected. 
testVolume = zeros(volumeSize);

testVolume(curveIndexList(mainInVolInds)) = 1;

tempCC = bwconncomp(testVolume, 26);

tempStats = regionprops(tempCC, 'PixelList');

nRegions = length(tempStats);

%%% Get regions for top and bottom to add warning.

% plot regions for now 
figure; hold on; axis equal

cols = lines(nRegions);

for iRegion = 1:nRegions
    plot3(tempStats(iRegion).PixelList(:,1), tempStats(iRegion).PixelList(:,2), ...
        tempStats(iRegion).PixelList(:,3), '.');
    
end

%clear smallGrainExterior, clear smallGrainVolume, 
%clear loopVolume, clear centreCurveVolume
%% Test plot.
figure; hold on; axis equal; set(gca, 'Clipping', 'off')

line(xTopOfLoop*[1 1]+ xBoundsNew(1) - 1, [1 yTopOfLoop]+ yBoundsNew(1) - 1,...
    zTopOfLoop*[1 1]+ zBoundsNew(1) - 1)

line(xBottomOfLoop*[1 1]+ xBoundsNew(1) - 1, [1 yBottomOfLoop]+ yBoundsNew(1) - 1,...
    zBottomOfLoop*[1 1]+ zBoundsNew(1) - 1)

plot3(curveSubscriptArray(mainCurveInds,1), curveSubscriptArray(mainCurveInds,2), ...
   curveSubscriptArray(mainCurveInds,3), 'b.')

plot3(curveSubscriptArray(topCurveInds,1), curveSubscriptArray(topCurveInds,2), ...
    curveSubscriptArray(topCurveInds,3), 'r.')

plot3(curveSubscriptArray(endCurveInds,1), curveSubscriptArray(endCurveInds,2), ...
    curveSubscriptArray(endCurveInds,3), 'k.')

plot3(curveSubscriptArray(mainInVolInds,1), curveSubscriptArray(mainInVolInds,2), ...
    curveSubscriptArray(mainInVolInds,3), 'mx')

plot3(loopSubscriptArray(:,1)+ xBoundsNew(1) - 1, loopSubscriptArray(:,2)+ yBoundsNew(1) - 1,...
   loopSubscriptArray(:,3)+ zBoundsNew(1) - 1, 'go'); 

figure; hold on; axis equal; set(gca, 'Clipping', 'off')

plot3(curveSubscriptArray(mainInVolInds,1), curveSubscriptArray(mainInVolInds,2), ...
    curveSubscriptArray(mainInVolInds,3), 'mx')
%% Remove cut up from both exterior and main volume

%Just remove from volume where it overlap exterior
grainVolumeAlignedCut = grainVolumeAligned;

tempIndsTop = find(grainExterior(curveIndexList(topCurveInds)) == 1);
grainVolumeAlignedCut(curveIndexList(topCurveInds(tempIndsTop))) = 0;

tempIndsEnd = find(grainExterior(curveIndexList(endCurveInds)) == 1);
grainVolumeAlignedCut(curveIndexList(endCurveInds(tempIndsEnd))) = 0;

curveCutVolume = zeros(volumeSize, 'logical');

curveCutVolume(curveIndexList(topCurveInds(tempIndsTop))) = 1;

curveCutVolume(curveIndexList(endCurveInds(tempIndsEnd))) = 1;

grainExterior(curveIndexList(topCurveInds(tempIndsTop))) = 0;

grainExterior(curveIndexList(endCurveInds(tempIndsEnd))) = 0;

figure; imshow(sum(grainExterior(:,:,1623),3))

%% Get surface of aleurone.
aleuroneExterior = grainExterior & (grainVolumeAligned == ALEURONE_INDEX);

% Get index list
aleuroneSurfaceIndexList = find(aleuroneExterior);

nIndex = length(aleuroneSurfaceIndexList); aleuroneSurfaceSubscriptArray = zeros(nIndex, 3);

[aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), aleuroneSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, aleuroneSurfaceIndexList);

%Get subscripts and volume of interior aleurone surface
% note, may not be fully connected

aleuroneInterior = grainVolumeAligned;

aleuroneInterior(aleuroneInterior == ALEURONE_INDEX) = 0;

% Aleruone interior label is currently endosperm and germ.
aleuroneInterior = logical(aleuroneInterior);

% Exterior is outer volume, and grow into other volume.
aleuroneInterior = imdilate(aleuroneInterior, STREL_18_CONNECTED);

% Previously check interior didn't overlap exterior, however, causes problems for thin aleurone where voxels may have both IDs
aleuroneInterior = (grainVolumeAligned == ALEURONE_INDEX) & aleuroneInterior;

% Get aleurone interior surface subscripts.
aleuroneInteriorIndexList = find(aleuroneInterior);

nIndex = length(aleuroneInteriorIndexList); aleuroneInteriorSubscriptArray = zeros(nIndex, 3);

[aleuroneInteriorSubscriptArray(:,1), aleuroneInteriorSubscriptArray(:,2), aleuroneInteriorSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, aleuroneInteriorIndexList);

%% Get germ surfaces
% Take largest region for germ.
%%% Shifted to using cut volume
germExterior = (grainExterior | curveCutVolume)  & (grainVolumeAligned == GERM_INDEX);

% Take largest connected region of surface.
tempCC = bwconncomp(germExterior, 26);

tempStats = regionprops(tempCC, 'PixelIdxList');

% Get number of voxels in each region. 
nRegions = length(tempStats); voxelsPerRegionArray = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
end

% Largest will generally be much larger than others.
[~, tempIndex] = max(voxelsPerRegionArray); tempStats(tempIndex) = [];

% Remove other regions from volume.
for iRegion = 1:nRegions-1
    germExterior(tempStats(iRegion).PixelIdxList) = 0;
end

germSurfaceIndexList = find(germExterior);

nIndex = length(germSurfaceIndexList); germSurfaceSubscriptArray = zeros(nIndex, 3);

[germSurfaceSubscriptArray(:,1), germSurfaceSubscriptArray(:,2), germSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, germSurfaceIndexList);

clear germExterior
%% Get exterior endosperm surface except for crease
endospermExterior = grainExterior & (grainVolumeAligned == ENDOSPERM_INDEX);

% Take tips of aleurone points.
curveCentre = mean(curveSubscriptArray(:,1));

tempInd = find(aleuroneSurfaceSubscriptArray(:,1) < curveCentre);

[~, topLeftTipInd] = max(aleuroneSurfaceSubscriptArray(tempInd,3));

topLeftTipInd = tempInd(topLeftTipInd);

tempInd = find(aleuroneSurfaceSubscriptArray(:,1) > curveCentre);

[~, topRightTipInd] = max(aleuroneSurfaceSubscriptArray(tempInd,3));

topRightTipInd = tempInd(topRightTipInd);

% Find closest points on aleurone surface.
[~, nearestGermLeft] = min(sqrt((germSurfaceSubscriptArray(:,1) - aleuroneSurfaceSubscriptArray(topLeftTipInd,1)).^2 + ...
    (germSurfaceSubscriptArray(:,2) - aleuroneSurfaceSubscriptArray(topLeftTipInd,2)).^2 + ...
    (germSurfaceSubscriptArray(:,3) - aleuroneSurfaceSubscriptArray(topLeftTipInd,3)).^2 )); 

[~, nearestGermRight] = min(sqrt((germSurfaceSubscriptArray(:,1) - aleuroneSurfaceSubscriptArray(topRightTipInd,1)).^2 + ...
    (germSurfaceSubscriptArray(:,2) - aleuroneSurfaceSubscriptArray(topRightTipInd,2)).^2 + ...
    (germSurfaceSubscriptArray(:,3) - aleuroneSurfaceSubscriptArray(topRightTipInd,3)).^2 ));

% Draw line between both sets of points, and remove from endosperm exterior

% Dilate line to ensure surface is cut.
dilateRadius = 2;

% Left side.
dMapFromGerm = bwdistgeodesic(grainExterior, sub2ind(volumeSize, ...
    aleuroneSurfaceSubscriptArray(topLeftTipInd,1), aleuroneSurfaceSubscriptArray(topLeftTipInd,2),...
    aleuroneSurfaceSubscriptArray(topLeftTipInd,3)), 'quasi-euclidean');

dMapFromAl = bwdistgeodesic(grainExterior, sub2ind(volumeSize, ...
    germSurfaceSubscriptArray(nearestGermLeft,1), germSurfaceSubscriptArray(nearestGermLeft,2),...
    germSurfaceSubscriptArray(nearestGermLeft,3)), 'quasi-euclidean');

dMap = dMapFromGerm + dMapFromAl;

dMap = round(dMap * 64)/64;

dMap(isnan(dMap)) = Inf;

lineVolume = imregionalmin(dMap);

lineVolume = imdilate(lineVolume, strel('sphere',dilateRadius)); 

leftLineInds = find(lineVolume);

nIndex = length(leftLineInds); leftLineSubscripts = zeros(nIndex, 3);

[leftLineSubscripts(:,1), leftLineSubscripts(:,2), leftLineSubscripts(:,3)] = ...
    ind2sub(volumeSize, leftLineInds);

% Right side.
dMapFromGerm = bwdistgeodesic(grainExterior, sub2ind(volumeSize, ...
    aleuroneSurfaceSubscriptArray(topRightTipInd,1), aleuroneSurfaceSubscriptArray(topRightTipInd,2),...
    aleuroneSurfaceSubscriptArray(topRightTipInd,3)), 'quasi-euclidean');

dMapFromAl = bwdistgeodesic(grainExterior, sub2ind(volumeSize, ...
    germSurfaceSubscriptArray(nearestGermRight,1), germSurfaceSubscriptArray(nearestGermRight,2),...
    germSurfaceSubscriptArray(nearestGermRight,3)), 'quasi-euclidean');

dMap = dMapFromGerm + dMapFromAl;

dMap = round(dMap * 64)/64;

dMap(isnan(dMap)) = Inf;

lineVolume = imregionalmin(dMap);

lineVolume = imdilate(lineVolume, strel('sphere',3));

rightLineInds = find(lineVolume);

nIndex = length(rightLineInds); rightLineSubscripts = zeros(nIndex, dilateRadius);

[rightLineSubscripts(:,1), rightLineSubscripts(:,2), rightLineSubscripts(:,3)] = ...
    ind2sub(volumeSize, rightLineInds);

clear dMapFromGerm, clear dMapFromAl, clear dMap, clear lineVolume;

% Cut from endosperm exterior.
endospermExterior(leftLineInds) = 0;

endospermExterior(rightLineInds) = 0;

% Split endosperm parts, remove largest that extends to bottom third of grain
tempCC = bwconncomp(endospermExterior, 26);

tempStats = regionprops(tempCC, 'PixelIdxList', 'PixelList');

% Get number of voxels in each region. 
nRegions = length(tempStats); voxelsPerRegionArray = zeros(nRegions,1);

minZPerRegion = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
    
    minZPerRegion(iRegion) = min(tempStats(iRegion).PixelList(:,3));
end

% Sort from large to small
[~, sortLength] = sort(voxelsPerRegionArray, 'descend');

% Sort by minimum Z value in 4 largest regions
[~, sortMinZ] = sort(minZPerRegion(sortLength(1:4)) ,'ascend');

% Remove 2 with lowest Z value but add to seperate volume
endospermCreaseExterior = zeros(volumeSize, 'logical');

figure; hold on; axis equal
for iRegion = 1:4
    plot3(tempStats(sortLength(sortMinZ(iRegion))).PixelList(:,1), tempStats(sortLength(sortMinZ(iRegion))).PixelList(:,2), ...
        tempStats(sortLength(sortMinZ(iRegion))).PixelList(:,3), '.')
end
title(sprintf('Will remove lowest %i', removeLowest));

%%% Would be good to add a test to see which large regions to consider, 
%%% and which low ones to remove

for iRegion = 1:removeLowest
    endospermExterior(tempStats(sortLength(sortMinZ(iRegion))).PixelIdxList) = 0;
    
    endospermCreaseExterior(tempStats(sortLength(sortMinZ(iRegion))).PixelIdxList) = 1;
end

%%% Note that some fluff along crease remains.
%%% Could grow on surface to connect these to main crease?

endospermSurfaceIndexList = find(endospermExterior);

nIndex = length(endospermSurfaceIndexList); endospermSurfaceSubscriptArray = zeros(nIndex, 3);

[endospermSurfaceSubscriptArray(:,1), endospermSurfaceSubscriptArray(:,2), endospermSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, endospermSurfaceIndexList);

% Get crease border to aleurone (not any endosperm or germ alonge this region)
endospermCreaseExteriorBorder = imdilate(endospermCreaseExterior, STREL_18_CONNECTED) & aleuroneExterior;

% plot
figure; hold on ; axis equal; set(gca, 'Clipping', 'off')

plot3(aleuroneSurfaceSubscriptArray(1:100:end,1), aleuroneSurfaceSubscriptArray(1:100:end,2), ...
     aleuroneSurfaceSubscriptArray(1:100:end,3), 'b.')
 
plot3(endospermSurfaceSubscriptArray(:,1), endospermSurfaceSubscriptArray(:,2), endospermSurfaceSubscriptArray(:,3), 'g.')

plot3(germSurfaceSubscriptArray(:,1), germSurfaceSubscriptArray(:,2), germSurfaceSubscriptArray(:,3), 'y.')

plot3(leftLineSubscripts(:,1), leftLineSubscripts(:,2), leftLineSubscripts(:,3), 'kx')

plot3(rightLineSubscripts(:,1), rightLineSubscripts(:,2), rightLineSubscripts(:,3), 'kx')

clear endospermCreaseExterior
%% Check that aleurone and endosperm exteriors form continous region
% Not required for either indvidually, but should work for whole

combinedExterior = endospermExterior + aleuroneExterior;

% Edge will be fully connected by defeault.
%combinedExterior(aleuroneEdgeIndexList) = 3;

% Take largest connected region.
tempCC = bwconncomp(combinedExterior, 26);

tempStats = regionprops(tempCC, 'PixelIdxList');

% Get number of voxels in each region. 
nRegions = length(tempStats); voxelsPerRegionArray = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
end

% Largest will generally be much larger than others.
[~, tempIndex] = max(voxelsPerRegionArray); 

tempStats(tempIndex) = []; voxelsPerRegionArray(tempIndex) = [];

% Remove other regions from volumes and collect indices to remove from lists.
endoToRemove = zeros(sum(voxelsPerRegionArray),1)*NaN;

aleuroneToRemove = zeros(sum(voxelsPerRegionArray),1)*NaN;

counter = 1;

for iRegion = 1:nRegions-1
    
    combinedExterior(tempStats(iRegion).PixelIdxList) = 0;
    
    for jIndex = 1:length(tempStats(iRegion).PixelIdxList)
        
        % Get values.
        tempIndex = tempStats(iRegion).PixelIdxList(jIndex);
        
        value = grainVolumeAligned(tempIndex);
        
        %Delete from appropriate volume and record for list
        if value == ENDOSPERM_INDEX
            %endosperm
            endospermExterior(tempIndex) = 0;
            
            endoToRemove(counter) = find(endospermSurfaceIndexList == tempIndex);
                        
        elseif value == ALEURONE_INDEX
            %aleurone surfaces
            aleuroneExterior(tempIndex) = 0;
            
            aleuroneToRemove(counter) = find(aleuroneSurfaceIndexList == tempIndex);
            
        end
        
        counter = counter + 1;
    end
end

%remove from lists and arrays
endoToRemove(isnan(endoToRemove)) = []; aleuroneToRemove(isnan(aleuroneToRemove)) = [];

endospermSurfaceIndexList(endoToRemove) = [];

endospermSurfaceSubscriptArray(endoToRemove,:) = [];

aleuroneSurfaceIndexList(aleuroneToRemove) = [];

aleuroneSurfaceSubscriptArray(aleuroneToRemove,:) = [];

%% Calculate edge of combined surface.
%%% Changed to cut volume
combinedEdge = imdilate((grainExterior | curveCutVolume) & ~combinedExterior, STREL_18_CONNECTED);

combinedEdge = combinedEdge & combinedExterior;

% Take largest connected region of edge.
tempCC = bwconncomp(combinedEdge, 26);

tempStats = regionprops(tempCC, 'PixelIdxList');

% Get number of voxels in each region. 
nRegions = length(tempStats); voxelsPerRegionArray = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
end

% Largest will generally be much larger than others.
[~, tempIndex] = max(voxelsPerRegionArray); tempStats(tempIndex) = [];

% Remove other regions from volume.
for iRegion = 1:nRegions-1
    combinedEdge(tempStats(iRegion).PixelIdxList) = 0;
end

combinedEdgeIndexList = find(combinedEdge);

nIndex = length(combinedEdgeIndexList); combinedEdgeSubscriptArray = zeros(nIndex, 3);

[combinedEdgeSubscriptArray(:,1), combinedEdgeSubscriptArray(:,2), combinedEdgeSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, combinedEdgeIndexList);

creaseAleuroneInds = find(endospermCreaseExteriorBorder(combinedEdgeIndexList));

clear combinedExterior, clear combinedEdge, %clear curveCutVolume

%% Calculate surface area of aleurone exterior and interior
%https://se.mathworks.com/matlabcentral/answers/93023-is-there-a-matlab-function-that-can-compute-the-area-of-my-patch
% Matched Amira results closely for bee eye, should test for grain

warning('Test area calc matches Amira')

%Exterior
tic
tempSurf = isosurface(aleuroneExterior,0.5);
toc
verts = tempSurf.vertices*VOXEL_SIZE;

faces = tempSurf.faces;

a = verts(faces(:, 2), :) - verts(faces(:, 1), :);

b = verts(faces(:, 3), :) - verts(faces(:, 1), :);

c = cross(a, b, 2);

aleuroneExteriorArea = 1/2 * sum(sqrt(sum(c.^2, 2)))/2;

%Interior
tic
tempSurf = isosurface(aleuroneInterior,0.5);
toc

verts = tempSurf.vertices*VOXEL_SIZE;

faces = tempSurf.faces;

a = verts(faces(:, 2), :) - verts(faces(:, 1), :);

b = verts(faces(:, 3), :) - verts(faces(:, 1), :);

c = cross(a, b, 2);

aleuroneInteriorArea = 1/2 * sum(sqrt(sum(c.^2, 2)))/2;

clear tempSurf

%% Allocate points equally across surface and border, as in bee
%Create combined surface of endosperm and aleurone

%Id is coded 1 for aleurone, 0 for endosperm
combinedSurfaceSubscripts = [aleuroneSurfaceSubscriptArray' endospermSurfaceSubscriptArray']';

combinedIdentity = zeros(size(combinedSurfaceSubscripts,1), 1);

combinedIdentity(1:size(aleuroneSurfaceSubscriptArray,1)) = 1;

edgeIdentity = aleuroneExterior(combinedEdgeIndexList);

%To control choosing points
edgePointsToChoose = 1:length(combinedEdgeIndexList);

surfacePointsToChoose = 1:size(combinedSurfaceSubscripts,1);

sparsePointsToChoose = 1:size(combinedSurfaceSubscripts,1);

edgePointsChoosen = zeros(length(combinedEdgeIndexList), 1, 'logical');

surfacePointsChoosen = zeros(size(combinedSurfaceSubscripts,1), 1, 'logical');

sparsePointsChoosen = zeros(size(combinedSurfaceSubscripts,1), 1, 'logical');

% On OB5
% For 20, 50, gives 343 1314 points 
% Seems good to get dense on border and moderate on surface

% Test edge points first. Select from top of germ (max Y)
tic
while ~isempty(edgePointsToChoose)
    [~, ind] = max( combinedEdgeSubscriptArray(edgePointsToChoose,3));
    
    pointChoosen = combinedEdgeSubscriptArray(edgePointsToChoose(ind),:);
    
    edgePointsChoosen(edgePointsToChoose(ind)) = 1;
    
    % Remove other edge points within edge distance
    distances = sqrt( (combinedEdgeSubscriptArray(edgePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (combinedEdgeSubscriptArray(edgePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (combinedEdgeSubscriptArray(edgePointsToChoose,3) - pointChoosen(3)).^2);
    
    edgePointsToChoose(distances < edgeDistance) = [];
    
    % Remove surface points within edge distance
    distances = sqrt( (combinedSurfaceSubscripts(surfacePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (combinedSurfaceSubscripts(surfacePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (combinedSurfaceSubscripts(surfacePointsToChoose,3) - pointChoosen(3)).^2);
    
    surfacePointsToChoose(distances < edgeDistance) = [];
end
toc

edgePointsChoosen = find(edgePointsChoosen); 

[~, creasePointsInEdge] = intersect(edgePointsChoosen, creaseAleuroneInds);

% Select surface points from remaining, again going down Z
tic
while ~isempty(surfacePointsToChoose)
    [~, ind] = max(combinedSurfaceSubscripts(surfacePointsToChoose,3));
    
    pointChoosen = combinedSurfaceSubscripts(surfacePointsToChoose(ind),:);
    
    surfacePointsChoosen(surfacePointsToChoose(ind)) = 1;
    
    % Remove other surface points within surface distance
    distances = sqrt( (combinedSurfaceSubscripts(surfacePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (combinedSurfaceSubscripts(surfacePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (combinedSurfaceSubscripts(surfacePointsToChoose,3) - pointChoosen(3)).^2);
    
    surfacePointsToChoose(distances < surfaceDistance) = [];
end
toc

surfacePointsChoosen = find(surfacePointsChoosen); 

% Select sparse points by z layer for speed

sparseLayers = min(combinedSurfaceSubscripts(:,3)):sparsePointsDistance: ...
    max(combinedSurfaceSubscripts(:,3));
tic
for iLayer = sparseLayers
    
    pointsOnLayer = find(combinedSurfaceSubscripts(sparsePointsToChoose,3) == iLayer);
    
    while ~isempty(pointsOnLayer)
        %Select from max Y down
        %[~, ind] = max(combinedSurfaceSubscripts(sparsePointsToChoose(pointsOnLayer),2));
        
        % Take random point on layer
        ind = round(rand(1)*(length(pointsOnLayer)-1)) + 1;
        
        pointChoosen = combinedSurfaceSubscripts(sparsePointsToChoose(pointsOnLayer(ind)),:);

       % Use different test scheme for sparse points to try and improve speed.
       % Check if point less than sparse distance to any previously choosen edge or surface points.
       % Note, moved to 2D.
        testEdge = sum( sqrt( (combinedEdgeSubscriptArray(edgePointsChoosen,1) - pointChoosen(1)).^2 + ...
            (combinedEdgeSubscriptArray(edgePointsChoosen,2) - pointChoosen(2)).^2) < sparsePointsDistance) == 0;

        testSurface = sum( sqrt( (combinedSurfaceSubscripts(surfacePointsChoosen,1) - pointChoosen(1)).^2 + ...
            (combinedSurfaceSubscripts(surfacePointsChoosen,2) - pointChoosen(2)).^2) < sparsePointsDistance) == 0;

        if testEdge && testSurface
            % If not, test distance to exisiting sparse points - seems inefficient

            sparsePointsChoosen(sparsePointsToChoose(pointsOnLayer(ind))) = 1;

            % Remove other surface points within surface distance from both lists
            % Note: now done in 2D
            indsOut = find(sqrt( (combinedSurfaceSubscripts(sparsePointsToChoose(pointsOnLayer),1) - pointChoosen(1)).^2 + ...
                (combinedSurfaceSubscripts(sparsePointsToChoose(pointsOnLayer),2) - pointChoosen(2)).^2) < sparsePointsDistance);
            
            pointsOnLayer(indsOut) = [];

        else
            %Remove point from both lists
            pointsOnLayer(ind) = [];
        end

        %Remove point
        %sparsePointsToChoose(ind) = [];

    end
end
toc

sparsePointsChoosen = find(sparsePointsChoosen);

[length(edgePointsChoosen) length(surfacePointsChoosen) length(sparsePointsChoosen)]

tesFig = figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
plot3(combinedSurfaceSubscripts(surfacePointsChoosen,1), combinedSurfaceSubscripts(surfacePointsChoosen,2), ...
    combinedSurfaceSubscripts(surfacePointsChoosen,3), 'gx')

plot3(combinedEdgeSubscriptArray(edgePointsChoosen,1), combinedEdgeSubscriptArray(edgePointsChoosen,2), ...
    combinedEdgeSubscriptArray(edgePointsChoosen,3), 'ro')

plot3(combinedSurfaceSubscripts(sparsePointsChoosen,1), combinedSurfaceSubscripts(sparsePointsChoosen,2), ...
    combinedSurfaceSubscripts(sparsePointsChoosen,3), 'm.')

plot3(endospermSurfaceSubscriptArray(1:1:end,1), endospermSurfaceSubscriptArray(1:1:end,2), ...
     endospermSurfaceSubscriptArray(1:1:end,3), 'b.')
%% Load grey image.
% Loading tiff stack takes a while.
greyVolume = loadtiffstack(greyDirectory, 1);

if any(padLow); greyVolume = padarray(greyVolume, padLow, 0, 'pre'); end

if any(padHigh); greyVolume = padarray(greyVolume, padHigh, 0, 'post'); end

if size(greyVolume) ~= volumeSize; error('Size mismatch'); end

greyVolumeAligned = uint8( affine_transform_full(single(greyVolume), M1*M2*M3*M4, 1));

% Affine transform with linear interp. does not seem to have small holes problem.

clear greyVolume
%% Calulate distance between map points

%%% Geodesic connectivity calculated with 26

subscriptsToInterpolate = [combinedSurfaceSubscripts(surfacePointsChoosen,:)'...
    combinedEdgeSubscriptArray(edgePointsChoosen,:)']';

indsToInterpolate = sub2ind(volumeSize, subscriptsToInterpolate(:,1), subscriptsToInterpolate(:,2),...
    subscriptsToInterpolate(:,3));

interpolatedIdentity = [combinedIdentity(surfacePointsChoosen)' edgeIdentity(edgePointsChoosen)']';

indsForSparse = sub2ind(volumeSize, combinedSurfaceSubscripts(sparsePointsChoosen,1), ...
    combinedSurfaceSubscripts(sparsePointsChoosen,2), combinedSurfaceSubscripts(sparsePointsChoosen,3));

subscriptsForSparse = combinedSurfaceSubscripts(sparsePointsChoosen,:);

sparseIdentity = combinedIdentity(sparsePointsChoosen);

nPoints = length(indsToInterpolate);

nSparsePoints = length(indsForSparse);

% Combine index lists for normal calculations.
totalIndexList = [aleuroneSurfaceIndexList' endospermSurfaceIndexList']';

% Make distance matrix single precision to save space (particularly on sparse and combined dMap)
distanceMatrix = zeros(nPoints, nPoints, 'single')*NaN;

distanceMatrixSparse = zeros(nPoints, nSparsePoints, 'single')*NaN;

indexsToUseByPoint = cell(nPoints,1);

% Set up for normal calculation
normalRadius = 50;

% Loop through geodesic distance calculations for each point then put into matrix.

parfor iPoint = 1:nPoints %
    tic
    % Try using grain exterior volume, can traverese germ, but not curve cut
    dMap = bwdistgeodesic(grainExterior, indsToInterpolate(iPoint),'quasi-euclidean');
    toc
    
    % Pause to let matlab free memory (?)
    pause(0.1)
    
    distanceMatrix(iPoint, :) = dMap(indsToInterpolate);
    
    distanceMatrixSparse(iPoint,:) = dMap(indsForSparse);
    
    % Select points based on distance map to prevent picking points on both sides of crease.
    indexsToUseByPoint{iPoint} = find(dMap(totalIndexList) < normalRadius*2);
end
  
save(sprintf('C:\\Users\\Admin\\Documents\\MATLAB\\Temp_data\\%s_temp', 'distanceMatrix'), ...
    'edgeDistance', 'surfaceDistance', 'sparsePointsDistance', 'normalRadius',...    
    'distanceMatrix', 'distanceMatrixSparse', 'indexsToUseByPoint', '-v7.3');

%% Calculate normals and thicknesses
% Find sparse points associated with each main point, these will be referenced for calculating normals.
%%% Note this has to be done using distances from map, if done using 3D
%%% distance points on opposite sides of crease can be linked
%%% So, required splitting for loop between distance map and normals

if nPoints ~= size(distanceMatrixSparse,1) | nSparsePoints ~= size(distanceMatrixSparse,2)
    error('Mismatched sizes and matrix - load wrong temp file')
end

pointToSparseLinks = cell(nPoints,1);

for iPoint = 1:nSparsePoints
    
    %Find closest main point
    [~, closestPointIndex] = min( distanceMatrixSparse(:, iPoint) );
    
    temp = pointToSparseLinks{closestPointIndex};
    
    temp = [temp iPoint];
    
    pointToSparseLinks{closestPointIndex} = temp;
end

normalByPoint = zeros(nPoints, 3)*NaN; normalForSparseCell = cell(nPoints, 1);  

internalIntersectByPoint = zeros(nPoints, 3)*NaN; internalIntersectForSparseCell = cell(nPoints, 1);

% Intenstiy and thickness may not be very useful for points on edge.

thicknessByPoint = zeros(nPoints, 1)*NaN; thicknessForSparseCell = cell(nPoints, 1); 

averageIntensityByPoint = zeros(nPoints, 1)*NaN; averageIntensityForSparseCell = cell(nPoints, 1);

grid3D.nx = volumeSize(1); grid3D.ny = volumeSize(2); grid3D.nz = volumeSize(3);
grid3D.minBound = [1 1 1]';
grid3D.maxBound = volumeSize';

% Calculate normals to get thickness
parfor iPoint = 1:nPoints %
%for iPoint = 25; %1761:nPoints
    % Make list with this location and sparse points
    sparseLinks = pointToSparseLinks{iPoint};
    
    % Get indexes to use
    indexListInRange = indexsToUseByPoint{iPoint};
    
    % Make temporary arrays for results.
    tempNormalArray = zeros(length(sparseLinks),3)*NaN;
    
    tempInternalIntersectArray = zeros(length(sparseLinks),3)*NaN;
    
    tempThicknessArray = zeros(length(sparseLinks),1)*NaN;
    
    tempAverageIntensityArray = zeros(length(sparseLinks),1)*NaN;
    
    subscriptsToCalculate = [subscriptsToInterpolate(iPoint,:)' ...
        subscriptsForSparse(sparseLinks,:)']';
    
    indsToCalculate = [indsToInterpolate(iPoint), indsForSparse(sparseLinks)']';
    
    pointIdentities = [interpolatedIdentity(iPoint) sparseIdentity(sparseLinks)']';
    
    for jSubscript = 1:size(subscriptsToCalculate,1) 
        currentSubscript = subscriptsToCalculate(jSubscript,:);
        
        % Check that target is included in patch to calculate
        if ~any(indsToCalculate(jSubscript) == totalIndexList(indexListInRange))
            error('%i %i - Target point not in patch', iPoint, jSubscript)
        end
            
        tempNormal = normaltosurface(currentSubscript , ...
            combinedSurfaceSubscripts(indexListInRange,:), [], [], normalRadius);
        
        % Test normal direction by looking in both directions for greater than 2 lengths of air
        depth = 10;
        
        [tempX, tempY, tempZ, forwardDistance] = amanatideswooalgorithm_efficient(currentSubscript, tempNormal, grid3D, 0,...
            [], depth, 1);
        
        % Tightest gap is likely to be on crease cut by curve (~1 pixels)
        forwardIndexList = sub2ind(volumeSize, tempX, tempY, tempZ);
         
        forwardSum = cumsum((grainVolumeAlignedCut(forwardIndexList) == 0 ) .* forwardDistance);
        
        forwardAirIndex = find(forwardSum > 0.5);
        
        [tempX, tempY, tempZ, backwardDistance] = amanatideswooalgorithm_efficient(currentSubscript, -tempNormal, grid3D, 0,...
            [], depth, 1);
        
        backwardIndexList = sub2ind(volumeSize, tempX, tempY, tempZ);
        
        backwardSum = cumsum((grainVolumeAlignedCut(backwardIndexList) == 0 ) .* backwardDistance);

        backwardAirIndex = find(backwardSum > 0.5);
        
        % [grainVolumeAlignedCut(forwardIndexList(1:depth)) grainVolumeAlignedCut(backwardIndexList(1:depth))]
        
        gotDirection = 0;
        
        % Pick direction based on indexes
        if isempty(forwardAirIndex) && ~isempty(backwardAirIndex)
            gotDirection = 1;
            
        elseif ~isempty(forwardAirIndex) & isempty(backwardAirIndex)
            gotDirection = 1;
            
            tempNormal = -tempNormal;
        elseif ~isempty(forwardAirIndex) & ~isempty(backwardAirIndex)
            %Find length of non-air inds after 1st 2 (at end of list, after 1st air point).
            
            %warning('%i %i - Test this point', iPoint, jSubscript)  
            
            forwardEndSum = sum((grainVolumeAlignedCut(forwardIndexList(forwardAirIndex(1):end)) > 0) .* ...
                forwardDistance(forwardAirIndex(1):end));
            
            backwardEndSum = sum((grainVolumeAlignedCut(backwardIndexList(backwardAirIndex(1):end)) > 0) .* ...
                backwardDistance(backwardAirIndex(1):end));
            
            % Correct if there is a clear direction
            % Note, this will mostly likely get discarded in intersection
            % and intergral test, but aim here is to get correct direction
            if forwardEndSum > 4*backwardEndSum
                gotDirection = 1;
                
            elseif 4*forwardEndSum < backwardEndSum
                 gotDirection = 1;
            
                tempNormal = -tempNormal;
            else
                % Cannot determine correct normal direction from this info
                
                warning('%i %i - No clear direction - some air inds on both sides', iPoint, jSubscript)

                gotDirection = 0;
            end
            %
            
        elseif isempty(forwardAirIndex) & isempty(backwardAirIndex)
            % Direction is not clear, and normal is not very meaningfull
            
            % Can occur if exerior is kind of indented and remainder of
                % surface pulls normal so it has not air intersect
            % Most likely to occur on edge points
            
            warning('%i %i - No clear direction - no air inds', iPoint, jSubscript)
            
            gotDirection = 0;
        end
        
        % Save normal if direction found
        if gotDirection
            if jSubscript == 1
                normalByPoint(iPoint,:) = tempNormal;
            else
                tempNormalArray(jSubscript-1,:) = tempNormal;
            end
        end
        
        % Continue if direction found and part of aleurone
        if gotDirection & pointIdentities(jSubscript)
 
            % Plot for testing - note line points out while normal actually points in.
            if testPlotNormals & jSubscript == 1
                figure(tesFig)
                if iPoint < surfacePointsChoosen
                    lineColor = 'g';
                else
                    lineColor = 'r';
                end
%                 figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
%                 plot3(combinedSurfaceSubscripts(surfacePointsChoosen,1), combinedSurfaceSubscripts(surfacePointsChoosen,2), ...
%                     combinedSurfaceSubscripts(surfacePointsChoosen,3), 'gx')
% % 
%                 plot3(combinedEdgeSubscriptArray(edgePointsChoosen,1), combinedEdgeSubscriptArray(edgePointsChoosen,2), ...
%                     combinedEdgeSubscriptArray(edgePointsChoosen,3), 'ro')
%                 
%                 plot3(combinedSurfaceSubscripts(indexListInRange,1), combinedSurfaceSubscripts(indexListInRange,2),...
%                     combinedSurfaceSubscripts(indexListInRange,3), '.')
% % 
%                  plot3(currentSubscript(1), currentSubscript(2), currentSubscript(2), 'c*')
% %                 
%                 line([0 tempNormal(1)*200]+currentSubscript(1),...
%                     [0 tempNormal(2)*200]+currentSubscript(2), ...
%                     [0 tempNormal(3)*200]+currentSubscript(3), 'color', 'r')

                
                figure; imshow(grainVolumeAlignedCut(:,:,currentSubscript(3))*100)
                hold on
                
%                 line( [0 tempNormal(2)*200]+currentSubscript(2),  ...
%                     [0 tempNormal(1)*200]+currentSubscript(1),...
%                         'color', 'r')
                    
                plot(currentSubscript(2), currentSubscript(1) , 'c*')

                plot(combinedSurfaceSubscripts(indexListInRange,2),...
                    combinedSurfaceSubscripts(indexListInRange,1), '.')
            end

            %Draw line in voxel space.
            [tempX, tempY, tempZ, voxelDistances] = amanatideswooalgorithm_efficient(currentSubscript, ...
                tempNormal, grid3D, 0, [], maxAleuroneThickness, 1);

            indexList = sub2ind(volumeSize, tempX, tempY, tempZ);
            
            % Test if start point can be shifted back for aleurone
            %%% Need to prepend lists otherwise shift in line will occur due to rounding in amanatides 
            % Does not allow any air gap
            [tempXNeg, tempYNeg, tempZNeg, tempDistancesNeg] = amanatideswooalgorithm_efficient(currentSubscript, -tempNormal, grid3D, 0,...
                [], 10, 1);

            tempIndexList = sub2ind(volumeSize, tempXNeg, tempYNeg, tempZNeg);

            % Take first non aleurone index
            notAleuroneID = find(grainVolumeAlignedCut(tempIndexList) ~= ALEURONE_INDEX);

           
            if ~isempty(notAleuroneID)
               % Move start point back if possible, usual is 2
               if notAleuroneID(1) > 2
                    tempX = [fliplr(tempXNeg(2:(notAleuroneID(1)-1))') tempX']';

                    tempY = [fliplr(tempYNeg(2:(notAleuroneID(1)-1))') tempY']';

                    tempZ = [fliplr(tempZNeg(2:(notAleuroneID(1)-1))') tempZ']';

                    indexList = [fliplr(tempIndexList(2:(notAleuroneID(1)-1))') indexList']';

                    voxelDistances = [fliplr(tempDistancesNeg(2:(notAleuroneID(1)-1))') voxelDistances']';
               end
            else
                error('%i %i - No precceding non-aleurone', iPoint, jSubscript)
            end
            
            % Find intersects
            % Using interior intersects some extra skimming problems
            lineIDs = grainVolumeAlignedCut(indexList);
            
            interiorPoints = find(aleuroneInterior(indexList));

            %[lineIDs(1:50) interiorList(1:50)]
            
            intersectsCleared = 0;
            
            % Test there is an interior intersect
            if isempty(interiorPoints)
                % If not, figure outwhat happened.
                    % Ray starts on al. exterior, can hit interior or leave by al. exterior
                    % If re-enters, must do so through exerior
                    
                % Note that line could pass through curve cut without hitting exterior, but seems unlikely     
                exteriorInds = find(grainExterior(indexList));
                
                exteriorIDs = lineIDs(exteriorInds);
                
                if any(exteriorIDs ~= ALEURONE_INDEX)
                   % Has hit other exteriors
                   %%% Could try to test if there is a just a small air gap above aleurone, but discard for now
                   interiorPoints = [];
                           
                   intersectsCleared = 1;

                   warning('%i %i - No interior intersects', iPoint, jSubscript)   
               else
                   % Only passed through aleurone
                   error('%i %i - Only hits aleurone exterior, problem with normal', iPoint, jSubscript) 
               end
            end
            
            % Check other things do not occur before first interect
            if ~isempty(interiorPoints)
                if any(lineIDs(1:interiorPoints(1)) ~= ALEURONE_INDEX)
                    % Can sometimes be small pockets of air, test as in following section
                    if all((lineIDs(1:interiorPoints(1)) == ALEURONE_INDEX) | ...
                            (lineIDs(1:interiorPoints(1)) == 0))
                        % Test total lenth of air pockets
                        testInds = find(lineIDs(1:interiorPoints(1)) == 0);
                        
                        % If greater than threshold, have to do something
                        if sum(voxelDistances(testInds)) > blockThickness/2
                           % Could change intersect, but arguable if it should be before or after air, so just remove.
                           interiorPoints = [];
                           
                           intersectsCleared = 1;
                           
                           warning('%i %i - Too much air before', iPoint, jSubscript) 
                        end    
                    else
                        % This rarely happens, but can occur if normal leaves via exterior, reenters via endosperm, 
                        % then hits interior intersect... wtf?
                        
                         interiorPoints = [];
                           
                         intersectsCleared = 1;
                        
                         warning('%i %i - Unexpected materials before intersect', iPoint, jSubscript)
                    end
                end
            end
            
            %Find breaks in interior intersects
            interiorDiff = diff(interiorPoints);
            
            interiorSteps = find(interiorDiff > 1);
            
            if ~isempty(interiorSteps)
                for kStep = 1:length(interiorSteps)
                    testStart = interiorPoints(interiorSteps(kStep))+1;
                    
                    testEnd = testStart - 1 + interiorDiff(interiorSteps(kStep))-1;
                    
                    testInds = testStart:testEnd;
                    
                    % Test if gap filled by aleruone
                    % Weird that gap between al interior border should be filled by al, 
                        % but occurs due to skimming.
                    if any(lineIDs(testInds) ~= ALEURONE_INDEX)
                        % Test what length is filled by other IDs
                        testInds = testInds(lineIDs(testInds) ~= ALEURONE_INDEX);
                        
                        if sum(voxelDistances(testInds)) > blockThickness/2
                            % Remove if more than half block of other points
                            interiorPoints(interiorSteps(kStep)+1:end) = [];
                            
                            break;
                        end
                    end
                end
            end
            
            if ~isempty(interiorPoints)
                % Points is ok, record and thickness + intensity.
                tempIntersect = [tempX(interiorPoints(end)), ...
                    tempY(interiorPoints(end)), tempZ(interiorPoints(end))]; 

                tempThickness = sum(voxelDistances(1:interiorPoints(end)))*VOXEL_SIZE;

                tempIntensity = wmean( double(greyVolumeAligned( indexList(1:interiorPoints(end)))),...
                    voxelDistances(1:interiorPoints(end)));

                if jSubscript == 1            
                    internalIntersectByPoint(iPoint,:) = tempIntersect;

                    thicknessByPoint(iPoint) = tempThickness;

                    averageIntensityByPoint(iPoint) = tempIntensity;

                else
                    tempInternalIntersectArray(jSubscript-1,:) = tempIntersect;

                    tempThicknessArray(jSubscript-1) = tempThickness;

                    tempAverageIntensityArray(jSubscript-1) = tempIntensity;

                end

                if testPlotNormals
                    plot3(tempIntersect(1), tempIntersect(2), ...
                        tempIntersect(3),'mo');
                end
                
            else
               if ~intersectsCleared
                    % Almost all cases of this resolved by allowing interior to overlap exterior on thin al 
                    warning('%i %i - No interior intersect', iPoint, jSubscript);
               end
            end 
        end
    end
    
    % Store sparse point info into cell
    normalForSparseCell{iPoint} = tempNormalArray;
    
    internalIntersectForSparseCell{iPoint} = tempInternalIntersectArray;
    
    thicknessForSparseCell{iPoint} = tempThicknessArray;
    
    averageIntensityForSparseCell{iPoint} = tempAverageIntensityArray;
end

% Unload data from cells into arrays
normalForSparse = zeros(nSparsePoints, 3);  

internalIntersectForSparse = zeros(nSparsePoints, 3);

thicknessForSparse = zeros(nSparsePoints, 1); 

averageIntensityForSparse = zeros(nSparsePoints, 1);

for iPoint = 1:nPoints
    sparseLinks = pointToSparseLinks{iPoint};
    
    if ~isempty(sparseLinks)
        normalForSparse(sparseLinks,:) = normalForSparseCell{iPoint};

        internalIntersectForSparse(sparseLinks,:) = internalIntersectForSparseCell{iPoint};

        thicknessForSparse(sparseLinks) = thicknessForSparseCell{iPoint};

        averageIntensityForSparse(sparseLinks) = averageIntensityForSparseCell{iPoint};
    end
end

save(sprintf('C:\\Users\\Admin\\Documents\\MATLAB\\Temp_data\\%s_%i_%i_%i_%i_wEndo_distOnFull_voxelFractions_updateSkimming_splitCalcs', dataName, ...
        edgeDistance, surfaceDistance, sparsePointsDistance, normalRadius), ...
    'edgeDistance', 'surfaceDistance', 'sparsePointsDistance', 'normalRadius',...    
    'distanceMatrix', 'subscriptsToInterpolate', 'interpolatedIdentity',... 
    'distanceMatrixSparse', 'subscriptsForSparse', 'sparseIdentity',...
    'normalByPoint', 'internalIntersectByPoint', 'thicknessByPoint', 'averageIntensityByPoint',...
    'normalForSparse', 'internalIntersectForSparse', 'thicknessForSparse', 'averageIntensityForSparse',...
    'pointToSparseLinks', 'indexsToUseByPoint', 'curveStep', 'nrbStep', '-v7.3');

%load('/Users/gavintaylor/Documents/Matlab/Temp_data/distanceMatrix_10_50_3_100.mat')

%% Test plot normals
%%% Mostly ok, except for borders and possibly crease

figure; hold on; axis equal;

for iPoint = 1:nPoints
    tempNormal = normalByPoint(iPoint,:);
    
    startPoint = subscriptsToInterpolate(iPoint,:);

    plot3(startPoint(1), startPoint(2), startPoint(3), 'x')    
       
    line([0 50]*tempNormal(1) + startPoint(1), [0 50]*tempNormal(2) + startPoint(2), ...
        [0 50]*tempNormal(3) + startPoint(3)); 
end

%% plot distance map to to test
figure; hold on; axis equal;

ind = find((subscriptsToInterpolate(:,1)) == 434 & (subscriptsToInterpolate(:,2) == 301) & ...
    (subscriptsToInterpolate(:,3) == 35));

if isempty(ind)
   [~, ind] = min(subscriptsToInterpolate(:,2)); 
end

fscatter3(subscriptsToInterpolate(:,1), subscriptsToInterpolate(:,2), ...
    subscriptsToInterpolate(:,3),  distanceMatrix(ind,:), jet(100))

plot3(subscriptsToInterpolate(ind,1), subscriptsToInterpolate(ind,2), ...
    subscriptsToInterpolate(ind,3), 'rx')

%% Calaculate integrals under surface.

pointIntensityProfile = zeros(nPoints, numberOfBlocks)*NaN;

sparseIntensityProfile = zeros(nSparsePoints, numberOfBlocks)*NaN;

pointProfileID = zeros(nPoints, numberOfBlocks,4)*NaN;

sparseProfileID = zeros(nSparsePoints, numberOfBlocks,4)*NaN;

% If flag is not set, will keep just ignore block w/ other values, not break on them
allowOverhangs = 1;

% Some catches included to correct for some points being lost if normal
% skims border of other structure. - but some have to be lost...

% Step through main points
for iPoint = 1:(nPoints + nSparsePoints)
    
    if iPoint <= nPoints
        tempNormal = normalByPoint(iPoint,:);
        
        pointIdentity = interpolatedIdentity(iPoint);
        
        %Check if aleurone or endosperm
        if pointIdentity
            %Aleurone, extend line in from interior surface
            startPoint = internalIntersectByPoint(iPoint,:);
        else
           % Endosperm, extend in from outer surface 
           startPoint = subscriptsToInterpolate(iPoint,:);
        end
    else
        tempNormal = normalForSparse(iPoint-nPoints,:);
        
        pointIdentity = sparseIdentity(iPoint-nPoints);
        
        % As above
        if pointIdentity
            %Aleurone, extend line in from interior surface
            startPoint = internalIntersectForSparse(iPoint-nPoints,:);
        else
           % Endosperm, extend in from outer surface 
           startPoint = subscriptsForSparse(iPoint-nPoints,:);
        end
    end
    
%                 figure; imshow(grainVolumeAlignedCut(:,:,startPoint(3))*100)
%                 hold on
%                 
%                 line( [0 tempNormal(2)*200]+startPoint(2),  ...
%                     [0 tempNormal(1)*200]+startPoint(1),...
%                         'color', 'r')
%                     
%                 plot(startPoint(2), startPoint(1) , 'c*')
    
    if ~any(isnan(tempNormal)) & ~any(isnan(startPoint))
        %Draw line in voxel space.
        [tempX, tempY, tempZ, voxelDistances] = amanatideswooalgorithm_efficient(startPoint, ...
            tempNormal, grid3D, 0, [], depthToCalculate*2, 1);

        indexList = sub2ind(volumeSize, tempX, tempY, tempZ);

        if ~pointIdentity
            % Add preceding point so endosperm starts one voxel back
            [tempXNeg, tempYNeg, tempZNeg, tempDistancesNeg] = amanatideswooalgorithm_efficient(startPoint, ...
                -tempNormal, grid3D, 0, [], 10, 1);

            tempIndexList = sub2ind(volumeSize, tempXNeg, tempYNeg, tempZNeg);

            % Take first not endosperm index
            notEndospermID = find(grainVolumeAlignedCut(tempIndexList) ~= ENDOSPERM_INDEX);

            if ~isempty(notEndospermID)
               % Move start point back if possible

                tempX = [fliplr(tempXNeg(2:(notEndospermID(1)))') tempX']';

                tempY = [fliplr(tempYNeg(2:(notEndospermID(1)))') tempY']';

                tempZ = [fliplr(tempZNeg(2:(notEndospermID(1)))') tempZ']';

                indexList = [fliplr(tempIndexList(2:(notEndospermID(1)))') indexList']';

                voxelDistances = [fliplr(tempDistancesNeg(2:(notEndospermID(1)))') voxelDistances']';
            else
                error('%i No precceding non-endosperm', iPoint)
            end
        end

        lineIDs = grainVolumeAlignedCut(indexList);

        % Record ID profile, just the most common in each block
        % Will do down entire profile depth, even bits not considered for intensity later
        
        lineIDsTemp = lineIDs(2:end);

        voxelDistancesTemp = (voxelDistances(2:end));
        
        voxelSum = cumsum(voxelDistancesTemp);
        
        tempIDProfile = zeros(numberOfBlocks,4)*NaN;
        
        % Mostly copied from intensity averaging section
        for jBlock = 1:numberOfBlocks
            blockSet = find(voxelSum < blockThickness & voxelDistancesTemp ~= 0);

            if ~isempty(blockSet)
                distanceSet = voxelDistancesTemp(blockSet);

                %Check if final index in block is incomplete and there are following points
                if voxelSum(blockSet(end)) < blockThickness & ...
                        blockSet(end) ~= length(lineIDsTemp)
                    % If so, add next index after final to block
                    blockSet(end+1) = blockSet(end)+1;

                    %Add approrpaite distance, and calcualte amount remaining
                    distanceSet(end+1) = blockThickness - voxelSum(blockSet(end-1));

                    finalRemaining = voxelSum(blockSet(end)) - blockThickness;
                else
                   %Block is perfect length
                   finalRemaining = 0; 
                end

                % Just include if thickness more than 1/2 of full blcok
                if sum(distanceSet) > blockThickness/2
                    % Find length of each material type in block
                    %%% Note harded coded number of materials and values

                    for kMat = 1:4
                        tempIDProfile(jBlock, kMat) = sum(...
                            distanceSet(lineIDsTemp(blockSet) == ...
                            (kMat - 1) ));
                    end
                end

                % Zero out remaining blocks and retake average
                voxelDistancesTemp(blockSet) = 0;

                if finalRemaining ~= 0
                    voxelDistancesTemp(blockSet(end)) = finalRemaining;
                end

                voxelSum = cumsum(voxelDistancesTemp);
            end
        end
        
        % Save into correct array
        if iPoint <= nPoints
            pointProfileID(iPoint, :, :) = tempIDProfile*VOXEL_SIZE;
        else
            sparseProfileID(iPoint - nPoints, :, :) = tempIDProfile*VOXEL_SIZE;
        end
        
        % Get again
        voxelSum = cumsum(voxelDistances);
        
        % Get endosperm voxels
        endospermPoints = find(lineIDs == ENDOSPERM_INDEX);

        if ~isempty(endospermPoints)

            % Endosperm should not be on first point
            if endospermPoints(1) == 1
               error('%i Endosperm point too early', iPoint) 
            end

            % If overhangs not allowed, will remove both endosperm with
            % things above and also inclusions of other labels
            if allowOverhangs      
                % Will not be adding points, so need to redo list.
                if ~isempty(endospermPoints)
                    endospermPoints = endospermPoints(1):endospermPoints(end);
                end
            else
                if pointIdentity
                    % First point should always be al.
                    if lineIDs(1) ~= ALEURONE_INDEX
                       error('%i First point is not al.', iPoint) 
                    end

                    % Check if non-al things occur before first endosperm
                    if any(lineIDs(2:endospermPoints(1)-1) ~= ALEURONE_INDEX)
                        % Often bits of air or small things beforehand

                        testInds = find(lineIDs(2:endospermPoints(1)-1) ~= ALEURONE_INDEX);

                        % All half block of other materials (following proper aleurone intersect)
                        if sum(voxelDistances(testInds)) > blockThickness/2
                            warning('%i other things between aleurone and endosperm', iPoint)

                            endospermPoints = [];
                        end
                    end
                end

            
                % Remove endopserm after first 'real' break
                endospermDiff = diff(endospermPoints);

                endospermSteps = find(endospermDiff > 1);

                % Fun note: This seems superficially simlar to switch debouncing

                % Normal note: By doing this correction for both aleurone thickness and
                    % endosperm intensity, points on border could be counted in both if
                    % 1 1 2 1 2 2 interlaced pattern. 
                    % However, should be avoided by taking aleurone first and then continuing from there

                if ~isempty(endospermSteps)
                    % Check if each step is less than half block size

                    toInsert = [];
                    for jStep = 1:length(endospermSteps)

                        testStart = endospermPoints(endospermSteps(jStep))+1;

                        testEnd = testStart - 1 + endospermDiff(endospermSteps(jStep))-1;

                        testInds = testStart:testEnd;

                        if sum(voxelDistances(testInds)) < blockThickness/2
                            %Remove if there is not a (nearly) full block after
                             testStart = endospermPoints(endospermSteps(jStep))+endospermDiff(endospermSteps(jStep))-1;

                             % Length should be slightly less than thickness
                             testInds = find(voxelSum > voxelSum(testStart) & ...
                                voxelSum < (voxelSum(testStart) + blockThickness) );

                             % Check if last ind needs adding
                             if sum(voxelDistances(testInds)) < ...
                                 blockThickness

                                 testInds(end+1) = testInds(end) + 1;
                             end

                             testInds(testInds > length(lineIDs)) = [];

                             % Allow up to one voxel of other material
                             %%% Previous test didn't allow any

                             testInds = testInds(lineIDs(testInds) ~= ENDOSPERM_INDEX);

                             %if all(lineIDs(testInds) == ENDOSPERM_INDEX)
                             if sum(voxelDistances(testInds)) < 1
                                 % Insert as endosperm points
                                 toInsert = [toInsert endospermPoints(endospermSteps(jStep))+1:...
                                     endospermPoints(endospermSteps(jStep)+1)-1];
                             else
                                 % Remove following points
                                 endospermPoints(endospermSteps(jStep)+1:end) = [];

                                 warning('%i Part removed on following test', iPoint)

                                 break;
                             end
                        else
                            % Remove without test
                            endospermPoints(endospermSteps(jStep)+1:end) = [];

                            warning('%i Part removed on length test', iPoint)

                            break
                        end
                    end

                    endospermPoints = sort([endospermPoints' toInsert]');
                end
            end

            % Initalize lists with endosperm points
            indexListTemp = indexList(endospermPoints);

            lineIDsTemp = lineIDs(endospermPoints);
            
            voxelDistancesTemp = (voxelDistances(endospermPoints));

            voxelSum = cumsum(voxelDistancesTemp);
            
            % Need at least one full block included
            if sum(voxelDistancesTemp) > blockThickness

                distanceIncluded = 0;

                tempProfile = zeros(1, numberOfBlocks)*NaN;

                % Step through blocks, averaging intensity
                for jBlock = 1:numberOfBlocks

                    % Take first set of blocks
                    blockSet = find(voxelSum < blockThickness & voxelDistancesTemp ~= 0);

                    if ~isempty(blockSet)
                        distanceSet = voxelDistancesTemp(blockSet);

                        %Check if final index in block is incomplete and there are following points
                        if voxelSum(blockSet(end)) < blockThickness & ...
                                blockSet(end) ~= length(indexListTemp)
                            % If so, add next index after final to block
                            blockSet(end+1) = blockSet(end)+1;

                            %Add approrpaite distance, and calcualte amount remaining
                            distanceSet(end+1) = blockThickness - voxelSum(blockSet(end-1));

                            finalRemaining = voxelSum(blockSet(end)) - blockThickness;
                        else
                           %Block is perfect length
                           finalRemaining = 0; 
                        end

                        if allowOverhangs
                            % A bit of hacky add on, 
                            % essentially, will NaN out any block with <50p endosperm
                                % These will be considered as 'other' blocks in interpolator
                            % aim was to match block overhangs result
                                % all of block averaged if <50p endo,
                                % however other blockoverhangs tests that less than 1 non-endo follows
                                % this will include regardless of following
                                
                            testInds = find(lineIDsTemp(blockSet) ~= ENDOSPERM_INDEX);
                            
                            if sum(distanceSet(testInds)) < blockThickness/2
                                controlValue = 1;
                            else
                                controlValue = NaN; %Block it
                            end
                        else
                           controlValue = 1; 
                        end
                        
                        % Just include if thickness more than 3/4 of full block
                        % My be less at end
                        if sum(distanceSet) > blockThickness*3/4;
                            if size(distanceSet,2) ~= 1
                               distanceSet = distanceSet'; 
                            end
                            
                            % Take average
                            tempProfile(jBlock) = wmean( ...
                                double(greyVolumeAligned(indexListTemp(blockSet)))*...
                                controlValue, distanceSet);

                            distanceIncluded = distanceIncluded + sum(distanceSet); 
                        end

                        % Zero out remaining blocks and retake sum
                        voxelDistancesTemp(blockSet) = 0;

                        if finalRemaining ~= 0
                            voxelDistancesTemp(blockSet(end)) = finalRemaining;
                        end

                        voxelSum = cumsum(voxelDistancesTemp);
                    end
                end

                % Save into correct array
                if iPoint <= nPoints
                    pointIntensityProfile(iPoint, :) = tempProfile;

                else
                    sparseIntensityProfile(iPoint - nPoints, :) = tempProfile;

                end

                if distanceIncluded < depthToCalculate
                   warning('%i Full depth not reached', iPoint); 
                end
            else
               warning('%i Less than one block', iPoint) 
            end
        else
           if ~any(lineIDs == GERM_INDEX)
                error('%i No endosperm or germ, wrong normal?', iPoint) 
           end
        end
    else
       warning('%i Either no normal or no start point, skipped completely', iPoint) 
    end
end

%% display demo slice
figure;

zSlice = 800;

imshow(greyVolumeAligned(:,:,zSlice)'); hold on

aleuroneExeteriorIn = find(aleuroneSurfaceSubscriptArray(:,3) == zSlice);

aleuroneInteriorIn = find(aleuroneInteriorSubscriptArray(:,3) == zSlice);

plot(aleuroneSurfaceSubscriptArray(aleuroneExeteriorIn,1), aleuroneSurfaceSubscriptArray(aleuroneExeteriorIn,2), 'y.')

plot(aleuroneInteriorSubscriptArray(aleuroneInteriorIn,1), aleuroneInteriorSubscriptArray(aleuroneInteriorIn,2), 'y.')

allPointsInSlice = find(subscriptsToInterpolate(:,3) >= zSlice - 3 & ...
    subscriptsToInterpolate(:,3) <= zSlice + 3);

sparsePointsInSlice = find(subscriptsForSparse(:,3) >= zSlice - 3 & ...
    subscriptsForSparse(:,3) <= zSlice + 3) + nPoints;

cols = winterCoolColours;

for iPoint = [allPointsInSlice' sparsePointsInSlice']
    
    if iPoint <= nPoints
        tempNormal = normalByPoint(iPoint,:);
        
        pointIdentity = interpolatedIdentity(iPoint);
        
        %Check if aleurone or endosperm
        if pointIdentity
            %Aleurone, extend line in from interior surface
            startPoint = internalIntersectByPoint(iPoint,:);
        else
           % Endosperm, extend in from outer surface 
           startPoint = subscriptsToInterpolate(iPoint,:);
        end
    else
        tempNormal = normalForSparse(iPoint-nPoints,:);
        
        pointIdentity = sparseIdentity(iPoint-nPoints);
        
        % As above
        if pointIdentity
            %Aleurone, extend line in from interior surface
            startPoint = internalIntersectForSparse(iPoint-nPoints,:);
        else
           % Endosperm, extend in from outer surface 
           startPoint = subscriptsForSparse(iPoint-nPoints,:);
        end
    end
    
    if ~any(isnan(tempNormal)) & ~any(isnan(startPoint))
        %Draw line in voxel space.
        [tempX, tempY, tempZ, voxelDistances] = amanatideswooalgorithm_efficient(startPoint, ...
            tempNormal, grid3D, 0, [], depthToCalculate*2, 1);

        indexList = sub2ind(volumeSize, tempX, tempY, tempZ);
        
        % Don't bother stepping back
        
        voxelSum = cumsum(voxelDistances);
        
         for jBlock = 1:numberOfBlocks
            blockSet = find(voxelSum < blockThickness*jBlock & voxelSum ~= 0);
            
            voxelSum(blockSet) = 0;
            
            plot(tempX(blockSet), tempY(blockSet), 'color', cols(jBlock,:))
         end
    end
end
%% Calculate unwrapping 

% Check for points with inf distance.
if any(isinf(distanceMatrix(:))) || any(isnan(distanceMatrix(:)))
   error('Can not have NaN or Inf in distance matrix')
end

% Enforce symmetry (should be ok...). Should also be squared
distanceMatrixTemp = ((distanceMatrix.^2 + (distanceMatrix.^2)')/2);

% Apparently CMD is equivlenet to PCA; but seems to give different results
[coeff, pointsUnwrapped] = pca(distanceMatrixTemp, 'Algorithm','eig',...
    'Centered','on','NumComponents',2);

surfaceIndexList = 1:length(surfacePointsChoosen);

edgeIndexList = (length(surfacePointsChoosen)+1):(length(surfacePointsChoosen)+size(edgePointsChoosen,1));

%Get full edge aleurone list
aleuroneEdgeIndexList = edgeIndexList( find( interpolatedIdentity(edgeIndexList)));

endospermEdgeIndexList = edgeIndexList( find( ~interpolatedIdentity(edgeIndexList)));

% Split into crease and other aleurone
%%% Should do this earlier and clear endospermCreaseExteriorBorder
aleuroneCreaseEdgeIndexList = aleuroneEdgeIndexList( find( ...
    endospermCreaseExteriorBorder(indsToInterpolate(aleuroneEdgeIndexList))));

aleuroneOtherEdgeIndexList = aleuroneEdgeIndexList(find( ...
    ~endospermCreaseExteriorBorder(indsToInterpolate(aleuroneEdgeIndexList))));

% Test plot.
% figure; axis equal; hold on
% plot3(subscriptsToInterpolate(endospermEdgeIndexList,1), subscriptsToInterpolate(endospermEdgeIndexList,2), ...
%     subscriptsToInterpolate(endospermEdgeIndexList,3), 'gx');
% 
% plot3(subscriptsToInterpolate(aleuroneOtherEdgeIndexList,1), subscriptsToInterpolate(aleuroneOtherEdgeIndexList,2), ...
%     subscriptsToInterpolate(aleuroneOtherEdgeIndexList,3), 'yx');
% 
% plot3(subscriptsToInterpolate(aleuroneCreaseEdgeIndexList,1), subscriptsToInterpolate(aleuroneCreaseEdgeIndexList,2), ...
%     subscriptsToInterpolate(aleuroneCreaseEdgeIndexList,3), 'kx');

% Do embedding by placing extra points into space
sparsePointsUnwrapped = zeros(size(subscriptsForSparse,1),2, 'single');

% Take column wise mean on orginal distance map 
coloumnMeansDistanceMatrix = mean(distanceMatrixTemp,2);

%embeddingPseudoinverse = pinv(pointsUnwrapped(:,1:2));

for iPoint = 1:size(subscriptsForSparse,1)
% From: https://stats.stackexchange.com/questions/368331/project-new-point-into-mds-space
%     sparsePointsUnwrapped(iPoint,:) = -0.5*pointsUnwrapped(:,1:2)\(distanceMatrixSparse(:,iPoint).^2 -...
%         coloumnMeansDistanceMatrix);
    
%     sparsePointsUnwrapped(iPoint,:) = -0.5*embeddingPseudoinverse'*(distanceMatrixSparse(:,iPoint).^2 -...
%         coloumnMeansDistanceMatrix);

    %This seems most reliable
    sparsePointsUnwrapped(iPoint,:) = coeff'*(distanceMatrixSparse(:,iPoint).^2 -...
         coloumnMeansDistanceMatrix);
end

% Adjust scale of sparse to fit full (for some reason coeffs are small)
% sparseHegiht = max(sparsePointsUnwrapped(:,2)) - min(sparsePointsUnwrapped(:,2));
% 
% unwrappedHeight = max(pointsUnwrapped(:,2)) - min(pointsUnwrapped(:,2));
% 
% sparsePointsUnwrapped = sparsePointsUnwrapped/sparseHegiht*unwrappedHeight;

% For scaling, get highest non-crease aleruone edge point in bottom third of volume
indexListBottom3rd = find( subscriptsToInterpolate(aleuroneOtherEdgeIndexList, 3) < volumeSize(3)/3 );

[~, indexHighest] = max( subscriptsToInterpolate(aleuroneOtherEdgeIndexList(indexListBottom3rd), 2) );

indexBelow = aleuroneOtherEdgeIndexList(indexListBottom3rd(indexHighest));

% Find non-crease aleurone and endosperm edge in top half of volume that is best aligned in X
indexAleuroneListTopHalf = find( subscriptsToInterpolate(aleuroneOtherEdgeIndexList, 3) > volumeSize(3)/2 );

indexEndospermListTopHalf = find( subscriptsToInterpolate(endospermEdgeIndexList, 3) > volumeSize(3)/2 );

[aleuroneDist, indexClosestAleurone] = min( abs( ...
    subscriptsToInterpolate(aleuroneOtherEdgeIndexList(indexAleuroneListTopHalf), 1) - ...
    subscriptsToInterpolate(indexBelow,1)));

[endospermDist, indexClosestEndosperm] = min( abs( ...
    subscriptsToInterpolate(endospermEdgeIndexList(indexEndospermListTopHalf), 1) - ...
    subscriptsToInterpolate(indexBelow,1)));

if aleuroneDist < endospermDist
   indexAbove = (aleuroneOtherEdgeIndexList(indexAleuroneListTopHalf(indexClosestAleurone)));
else
   indexAbove = (endospermEdgeIndexList(indexEndospermListTopHalf(indexClosestEndosperm)));
end

%Take current unwrapped distance for scaling
unwrappedDistance = sqrt((pointsUnwrapped((indexAbove),1) - pointsUnwrapped((indexBelow),1)).^2 +...
    (pointsUnwrapped((indexAbove),2) - pointsUnwrapped((indexBelow),2)).^2);

% Find closest sparse points on distance map and scale to match distances
[~, closestBelowInd] = min(distanceMatrixSparse(indexBelow,:));

[~, closestAboveInd] = min(distanceMatrixSparse(indexAbove,:));

sparseUnwrappedDistance = sqrt((sparsePointsUnwrapped((closestAboveInd),1) - sparsePointsUnwrapped((closestBelowInd),1)).^2 +...
    (sparsePointsUnwrapped((closestAboveInd),2) - sparsePointsUnwrapped((closestBelowInd),2)).^2);

%Scale sparse based on 3D distance (not distance map, as we don't have dist between sparse points)
distOnFull = sqrt((subscriptsToInterpolate((indexAbove),1) - subscriptsToInterpolate((indexBelow),1)).^2 +...
    (subscriptsToInterpolate((indexAbove),2) - subscriptsToInterpolate((indexBelow),2)).^2 + ...
    (subscriptsToInterpolate((indexAbove),3) - subscriptsToInterpolate((indexBelow),3)).^2);

distOnSparse = sqrt((subscriptsForSparse((closestAboveInd),1) - subscriptsForSparse((closestBelowInd),1)).^2 +...
    (subscriptsForSparse((closestAboveInd),2) - subscriptsForSparse((closestBelowInd),2)).^2 + ...
    (subscriptsForSparse((closestAboveInd),3) - subscriptsForSparse((closestBelowInd),3)).^2);

if distOnSparse/distOnFull < 0.99
    error('Matching sparse points may be incorrect')
end

% Center both along above below line
sparsePointsUnwrapped(:,2) = sparsePointsUnwrapped(:,2) - mean(sparsePointsUnwrapped([closestAboveInd closestBelowInd], 2));

sparsePointsUnwrapped(:,1) = sparsePointsUnwrapped(:,1) - mean(sparsePointsUnwrapped([closestAboveInd closestBelowInd], 1));

pointsUnwrapped(:,2) = pointsUnwrapped(:,2) - mean(pointsUnwrapped([indexAbove indexBelow], 2));

pointsUnwrapped(:,1) = pointsUnwrapped(:,1) - mean(pointsUnwrapped([indexAbove indexBelow], 1));

% Scale based on relative distances for surf and unwrapping
distanceOnSurf = distanceMatrix((indexAbove), (indexBelow));

sparsePointsUnwrapped = (sparsePointsUnwrapped/sparseUnwrappedDistance)*distanceOnSurf*(distOnSparse/distOnFull);

pointsUnwrapped = (pointsUnwrapped/unwrappedDistance)*distanceOnSurf;

% % Match angle.
unwrappedAngle = atan2(pointsUnwrapped((indexAbove),2) - pointsUnwrapped((indexBelow),2), ...
    pointsUnwrapped((indexAbove),1) - pointsUnwrapped((indexBelow),1));

% Which one is below zero may depend on embedding...
pointsUnwrapped = [pointsUnwrapped ones(size(pointsUnwrapped,1),1) ] * ...
    make_transformation_matrix([0 0], -unwrappedAngle-pi/2+pi);

pointsUnwrapped = pointsUnwrapped(:,1:2);

sparsePointsUnwrapped = [sparsePointsUnwrapped ones(size(sparsePointsUnwrapped,1),1) ] * ...
    make_transformation_matrix([0 0], -unwrappedAngle-pi/2+pi);

sparsePointsUnwrapped = sparsePointsUnwrapped(:,1:2);

figure; 
%subplot(1,2,1);
hold on; axis equal

plot(sparsePointsUnwrapped(:,1), sparsePointsUnwrapped(:,2), 'b.');

plot(pointsUnwrapped(surfaceIndexList,1), pointsUnwrapped(surfaceIndexList,2), 'kx');

plot(pointsUnwrapped(edgeIndexList,1), pointsUnwrapped(edgeIndexList,2), 'ro');

plot(pointsUnwrapped(edgeIndexList(creasePointsInEdge),1), ...
    pointsUnwrapped(edgeIndexList(creasePointsInEdge),2), 'm*'); 

plot(pointsUnwrapped((indexAbove),1), pointsUnwrapped((indexAbove),2), 'gd');

plot(pointsUnwrapped((indexBelow),1), pointsUnwrapped((indexBelow),2), 'cd');

%% Create nice image plotting

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

thicknessImage = (thicknessImage-0)/(80-0);

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

intensityImage = (intensityImage-100)/(220-100);

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

figure; 

colsID = [[0 0 0]; [1 0 1]; [1 1 0]; [0 1 0]; [0 0 1]];

for iBlock = 1:numberOfBlocks    
    subplot(numPlotRows,5,iBlock);
    
    tempMap = zeros(xRange , yRange, 'uint8');
    
    tempMap(inMap) = IDProfileInterp(:, iBlock)+1;
    
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

figure;

% Put into arrays and plot
for iBlock = 1:numberOfBlocks
   
    tempImage = zeros(xRange , yRange);
    
    % Have to solve w/ color image because of scaling
    tempImageCol = zeros(xRange , yRange, 3);
    
    indsToUse = find(IDProfileInterp(:, iBlock) == ENDOSPERM_INDEX);
    
    % Get points and wrap to range
    tempImage(inMap(indsToUse)) = (intensityProfileInterp(indsToUse, iBlock)-40)/(120-40);
    
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
    
    subplot(numPlotRows,5,iBlock)

    imshow( permute(tempImageCol, [2 1 3])); hold on;
    
    plot(endoBorderX, endoBorderY, 'g.')
    plot(aleuroneBorderX, aleuroneBorderY, 'y.')
    
    title(sprintf('%i - %i um', [(iBlock-1) iBlock]*blockThickness*VOXEL_SIZE));
end
    
figure
colormap(gray(100))
hcb = colorbar; set(hcb,'Ticks', [0 0.5 1], 'TickLabels', {'40','80','120'})
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
for iBlock = 1:numberOfBlocks
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)) & ~isnan(thicknessByPoint));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock)) & ~isnan(thicknessForSparse));
    
    alIntensity = [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']';
    
    blockIntensity = [pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)']';
    
    subplot(numPlotRows,5,iBlock)
    plot(alIntensity, blockIntensity, '.');
    
    rValue = corr(alIntensity, blockIntensity);
    title(sprintf('%i - %i, r = %.2f', [(iBlock-1) iBlock]*blockThickness*VOXEL_SIZE, rValue));
end

subplot(numPlotRows,5,6)
xlabel('Aleurone intensity (AU)')
ylabel('Endosperm intensity (AU)')

figure;

for iBlock = 1:numberOfBlocks
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)) & ~isnan(thicknessByPoint));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock)) & ~isnan(thicknessForSparse));
    
    alThickness = [thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']';
    
    blockIntensity = [pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)']';
    
    subplot(numPlotRows,5,iBlock)
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
xlim([0 60]); ylim([0 10])

subplot(1,3,2)
plot(verticalIntensityMean, (verticalCoords+vertOffset)*VOXEL_SIZE/1000);
xlabel('Average aleurone intensity (AU)')
title('Vertical profiles')
xlim([0 200]); ylim([0 10])

subplot(1,3,3)
plot(verticalCount*areaPerPoint, (verticalCoords+vertOffset)*VOXEL_SIZE/1000);
xlabel('Average aleurone area (mm2)')
xlim([0 10]); ylim([0 10])

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

figure;
for iBlock = 1:numberOfBlocks
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock)));
    
    % Calc certical, using bins from al
    [verticalIntensityMean] = calculateProfile([pointIntensityProfile(valuesToUse, iBlock)' sparseIntensityProfile(valuesToSparse, iBlock)']',...
        -[offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']', vertBins, []);
    
    subplot(1,2,1); hold on
    plot(verticalIntensityMean, (verticalCoords+vertOffset)*VOXEL_SIZE/1000, 'color', cols(iBlock,:));

    [horizontalIntensityMean] = calculateProfile([pointIntensityProfile(valuesToUse, iBlock)' sparseIntensityProfile(valuesToSparse, iBlock)']',...
        [offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']', hoirzBins, xCent);

    subplot(1,2,2); hold on
    plot(horizontalCoords*VOXEL_SIZE/1000, horizontalIntensityMean, 'color', cols(iBlock,:));
end

subplot(1,2,1); hold on
xlabel('Average endosperm intensity (AU)')
ylabel('Distance from base (mm)')
xlim([0 100]); ylim([0 10])
title('Vertical profile of endosperm layers')

subplot(1,2,2); hold on
ylabel('Average endosperm intensity (AU)')
xlabel('Distance from centreline (mm)')
ylim([0 100]); xlim([-5 5])
title('Horizontal profile of endosperm layers')
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
