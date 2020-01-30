% Test load and plot of grain stack.
clc; clear; close all

% Other: Om_1_7_test
targetDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Grain LU/Data/OB6_test/Labels';

% Loading tiff stack takes a while.
grainVolume = loadtiffstack(targetDirectory, 1);

%grainVolume = grainVolume(:,1:500,:);

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

ALEURONE_INDEX = 1; % Material index.
ENDOSPERM_INDEX = 2;
GERM_INDEX = 3;
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

transform2Up = matrix2rotatevectors([0, 1, 0], grainCreaseAxis*transform2Vertical);



% Apply transformation to grain in volume.
temp = zeros(4,4); temp(4,4) = 1;

M1 = make_transformation_matrix(grainCenter - volumeSize/2);

M2 = temp; M2(1:3,1:3) = transform2Vertical;

M3 = temp; M3(1:3,1:3) = transform2Up;

M4 = make_transformation_matrix(volumeSize/2 - grainCenter);

grainVolumeAligned = uint8(affine_transform_full(single(grainVolume), M1*M2*M3*M4, 5));

grainVolumeAligned = imclose(grainVolumeAligned, STREL_6_CONNECTED);

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
nRegions = length(tempStats);

voxelsPerRegionArray = zeros(nRegions,1);

for iRegion = 1:nRegions
    voxelsPerRegionArray(iRegion) = length(tempStats(iRegion).PixelIdxList);
end

% Largest will generally be much larger than others.
[~, tempIndex] = max(voxelsPerRegionArray);

tempStats(tempIndex) = [];

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
aleuroneSurfaceIndexList = find(grainExterior & grainVolume == ALEURONE_INDEX);

nIndex = length(aleuroneSurfaceIndexList); aleuroneSurfaceSubscriptArray = zeros(nIndex, 3);

[aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), aleuroneSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, aleuroneSurfaceIndexList);

clear grainIndexList, clear grainVolume

%figure; hold on; axis equal; set(gca, 'Clipping', 'off')

%plot3(grainSurfaceSubscriptArray(:,1), grainSurfaceSubscriptArray(:,2), grainSurfaceSubscriptArray(:,3), 'b.')
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
        
        if sum(grainVolumeAligned(jColumn,:,iSlice))
            % Find highest point.
            temp = find(grainVolumeAligned(jColumn,:,iSlice) == ALEURONE_INDEX | ...
                grainVolumeAligned(jColumn,:,iSlice) == ENDOSPERM_INDEX);
            
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



% Calc distance map from plane underneath grain.
%basePlaneVolume = zeros(smallVolumeSize, 'logical');

%basePlaneVolume(:,1,:) = 1;

% Make mask image for grain.
grainMask = logical(smallGrainVolume);

% Step through and cut above centre line;
for iSlice = 1:smallVolumeSize(3)
    
    if ~isnan( centreLineHeightProfile(iSlice+zBoundsNew(1)-1))
        
        cutAbove = centreLineHeightProfile(iSlice+zBoundsNew(1)-1)-yBoundsNew(1)+1+50;
        
        if cutAbove < smallVolumeSize(2)
            
            grainMask(:, cutAbove:end, iSlice) = 1;
        end 
    else
        % Cut above centre of mass of pixels.
        sliceIndexList = find(smallGrainVolume(:,:,iSlice));

        [~, ySubscriptList, ] = ind2sub(smallVolumeSize, sliceIndexList); 
        
        if ~isempty(ySubscriptList)
            cutAbove = round( mean(ySubscriptList));

            grainMask(:, cutAbove:end, iSlice) = 1;
        end
    end
end

% Calculate distance maps - note polarity of mask is opposite when using bwdist
%distanceFromPlane = bwdistgeodesic(~grainMask, basePlaneVolume, 'quasi-euclidean');
tic
distanceFromGrain = bwdist(grainMask, 'quasi-euclidean');
toc

% Need to set interior to max to prevent traversal. 
distanceFromGrain(grainMask) = Inf;

% Graydist takes a long time. Squared to force closer to grain.
tic
distanceFromTop = graydist(distanceFromGrain.^2, sub2ind(smallVolumeSize, ...
    xTopOfLoop, yTopOfLoop, zTopOfLoop), 'quasi-euclidean');
toc

tic
distanceFromBottom = graydist(distanceFromGrain.^2, sub2ind(smallVolumeSize, ...
    xBottomOfLoop, yBottomOfLoop, zBottomOfLoop),'quasi-euclidean');
toc

dMap = distanceFromTop + distanceFromBottom;

clear distanceFromTop, clear distanceFromBottom

% Round to lower precision to prevent floating point errors.
%%% Note, changed from 8 to 32 as loop doesn't reach end,
%%% but higher causes some external points to be included. 
%%% 32 closes but allows some floaters
dMap = round(dMap * 2^5)/2^5;

dMap(isnan(dMap)) = Inf;



% Reginal minima should define connection.
loopVolume = imregionalmin(dMap);

clear dMap

%% Take top of connected curve.
%%% Don't need to do clean above.

coordinatesOfLoop = zeros(smallVolumeSize(3)*2,3)*NaN;

coordinatesOfLoop(1,:) = [xTopOfLoop, yTopOfLoop, zTopOfLoop];

[testX, testY, testZ] = meshgrid(-1:1, -1:1, -1:1);

testX(14) = []; testY(14) = []; testZ(14) = [];

roughXCentre = (xTopOfLoop+xBottomOfLoop)/2;

distanceFromEnd = bwdistgeodesic(loopVolume, sub2ind(smallVolumeSize, ...
  xBottomOfLoop, yBottomOfLoop, zBottomOfLoop), 'quasi-euclidean');

% Step along loop towards end.
for iStep = 2:smallVolumeSize(3)*2
    % Get all distances to the adjacent points of previous point.
    distanceArray = zeros(26,1);
    
    for jPoint = 1:26
        
        distanceArray(jPoint) = distanceFromEnd(coordinatesOfLoop(iStep-1,1)+testX(jPoint), ...
               coordinatesOfLoop(iStep-1,2)+testY(jPoint), coordinatesOfLoop(iStep-1,3)+testZ(jPoint)); 
           
    end

    % If just one point, take shortest.
    tempIndex = find(distanceArray == min(distanceArray));
    
    if length(tempIndex) == 1
        
        coordinatesOfLoop(iStep,:) = coordinatesOfLoop(iStep-1,:) + ...
            [testX(tempIndex) testY(tempIndex) testZ(tempIndex)];
        
        %plot3(coordinatesOfLoop(iStep,1), coordinatesOfLoop(iStep,2), coordinatesOfLoop(iStep,3), 'kx');
        
    else
        % If more than 1 point, take highest.
        tempIndex2 = find(testY(tempIndex) == max(testY(tempIndex)));

        if length(tempIndex2) == 1
            
            coordinatesOfLoop(iStep,:) = coordinatesOfLoop(iStep-1,:) + ...
                [testX(tempIndex(tempIndex2)) testY(tempIndex(tempIndex2)) testZ(tempIndex(tempIndex2))];

            %plot3(coordinatesOfLoop(iStep,1), coordinatesOfLoop(iStep,2), coordinatesOfLoop(iStep,3), 'go');
            
        else
           %%% To add, if multiple high points, take closest to centreline 
           
           error('Multiple high points.');
           
        end 
    end
    
    % If end point reached, break early.
    if all((coordinatesOfLoop(iStep,:)-[xBottomOfLoop, yBottomOfLoop, zBottomOfLoop]) == 0)
        
        break;
        
    end
end

coordinatesOfLoop(isnan(coordinatesOfLoop(:,1)),:) = [];

clear distanceThroughLoop;

%Clear loop and add replace with thin values.
loopVolume(:) = 0;

tempIndexList = sub2ind(smallVolumeSize, coordinatesOfLoop(:,1), coordinatesOfLoop(:,2), coordinatesOfLoop(:,3));

loopVolume(tempIndexList) = 1;

loopIndexList = find(loopVolume);

nIndex = length(loopIndexList); loopSubscriptArray = zeros(nIndex, 3);
%% Now find centre-curve by slice.

centreCurveVolume = zeros(smallVolumeSize, 'logical');

% Transform distance from grain so centrelines are emphasized.
distanceNearGrain = distanceFromGrain;

distanceNearGrain(grainMask) = 0;

distanceNearGrain = -log(distanceNearGrain);

distanceNearGrain = distanceNearGrain - min(distanceNearGrain(:));

distanceNearGrain(grainMask) = Inf;

basePlane = zeros(smallVolumeSize(1:2), 'logical');

basePlane(:,1) = 1;
% Note that using single point frequently causes multiple lines to be drawn.
%basePlane(round(smallVolumeSize(1)/2),1) = 1;

missingSlices = zeros(smallVolumeSize(3),1);

for iSlice = 1:smallVolumeSize(3)
   
    % Check if exisiting crease points of crease on slice.
    indexList = find(loopVolume(:,:,iSlice));
    
    if ~isempty(indexList)
        
        loopSlice = loopVolume(:,:,iSlice);
        
        % Calculate distance maps to find minimum path from top of crease to base plane.
        distanceFromCrease = graydist(distanceNearGrain(:,:,iSlice), loopSlice, 'quasi-euclidean');
        
        distanceFromBase = graydist(distanceNearGrain(:,:,iSlice), basePlane, 'quasi-euclidean');
        
        dMap = distanceFromCrease + distanceFromBase;
        
        %%% Need to be careful setting precision to balance floaters with holes.
        %%% Lower value reduces floaters but can cause slices not to be connected at top or bottom. 
        %%% Using full base plane or point also influences best precision.
        dMap = round(dMap * 2^5)/2^5;
        
        dMap(isnan(dMap)) = Inf;
        
        % Regional minima is shortest path
        curveSlice = imregionalmin(dMap);

        % If regions are not linked, everything will be Inf
        curveSlice(curveSlice & isinf(dMap)) = 0;
        
        if iSlice == 21
            %temp = curveSlice + grainMask(:,:,iSlice);
            %figure; imshow(temp)
        end

        % Test that line links top and bottom.
        planeInCurve = sum(curveSlice(:) & basePlane(:));
        
        loopInCurve = sum(curveSlice(:) & loopSlice(:));

        % Using plane helps, but occasionally double connections are made. 
        % Test for this and they will be interpolated if occuring.
        
        if planeInCurve == 1 & loopInCurve
            % Floaters and breaks can still occur (!)
            % Check top and bottom are in main regions.
            tempCC = bwconncomp(curveSlice, 8);

            tempStats = regionprops(tempCC, 'PixelIdxList');

            nRegions = length(tempStats);
            
            % Remove region from curve slice if top/bottom not included
            for jRegion = 1:nRegions
                if sum(loopSlice(tempStats(jRegion).PixelIdxList)) == 0 && ...
                    sum(basePlane(tempStats(jRegion).PixelIdxList)) == 0
                    
                    curveSlice(tempStats(jRegion).PixelIdxList) = 0; 
                    
                    %figure; imshow(curveSlice);
                end
            end
            
            % Curves can sometimes be thick, so now using thinning
            % Seems more noticable when single point on base plane used.
            curveSlice = bwmorph(curveSlice, 'thin', inf);
            
            %Add slice of curve and loop to curve volume
            centreCurveVolume(:,:,iSlice) = curveSlice | loopSlice;

            curveIndexList = find(curveSlice);

            nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

            [curveSubscriptArray(:,1), curveSubscriptArray(:,2)] = ...
               ind2sub(smallVolumeSize(1:2), curveIndexList);
           
           % Even with precision correction it's possible holes exist
           if length(unique(curveSubscriptArray(:,2))) ~= max(curveSubscriptArray(:,2))
               % Could also be that real holes in X exist, but wont be
               % flagged here.
               
               % If not at least one point per slice, interpolate this slice.
               %%% Note any found points will not be moved. 
               missingSlices(iSlice) = 1;
               %figure; plot(curveSubscriptArray(:,1), curveSubscriptArray(:,2),'.')
           end    
        else
           %Slice cannot be calculated, just add exisiting loop 
           centreCurveVolume(:,:,iSlice) = loopSlice;
           
           missingSlices(iSlice) = 1;
        end
    end
end

% basePlaneOnSlice = zeros(smallVolumeSize(3),1);
% for iSlice = 1:smallVolumeSize(3)
%     basePlaneOnSlice(iSlice) = sum(centreCurveVolume(:,1,iSlice));
% end
% figure; plot(basePlaneOnSlice)

%% Interpolate missing slices.
missingSliceIndexList = find(missingSlices);

nMissingSlices = length(missingSliceIndexList);

% Get points to interpolate.
yToInterpolate = zeros(nMissingSlices, smallVolumeSize(2))*NaN;

zToInterpolate = zeros(nMissingSlices, smallVolumeSize(2))*NaN;

for iSlice = 1:nMissingSlices
    % Find exisiting points on curve (will include loop).
    tempIndexList = find(centreCurveVolume(:,:,missingSliceIndexList(iSlice)));

    nIndex = length(tempIndexList); tempSubscriptArray = zeros(nIndex, 2);

    [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(smallVolumeSize(1:2), tempIndexList);
    
    % Only interpolate up to highest Y on loop, and remove exisiting.
    temp = 1:max(tempSubscriptArray(:,2));
    
    temp(unique(tempSubscriptArray(:,2))) = [];

    yToInterpolate(iSlice,temp) = temp;
    
    zToInterpolate(iSlice,:) = missingSliceIndexList(iSlice);
    
end

% Remove excess interpolation points.
zToInterpolate(isnan(yToInterpolate)) = [];

yToInterpolate(isnan(yToInterpolate)) = [];

curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, curveIndexList);

xInterpolated = round(griddata(curveSubscriptArray(:,2), curveSubscriptArray(:,3), curveSubscriptArray(:,1), ...
    yToInterpolate, zToInterpolate, 'linear'));

% Add interpolated points to curve
tempIndex = sub2ind(smallVolumeSize, xInterpolated, yToInterpolate, zToInterpolate);

% Delete this and fix problem.
%centreCurveVolume(tempIndex(~isnan(tempIndex))) = 1;

% Test plot.
figure; hold on; axis equal; set(gca, 'Clipping', 'off')

%plot3(grainSurfaceSubscriptArray(:,1)-xBoundsNew(1)+1, grainSurfaceSubscriptArray(:,2)-yBoundsNew(1)+1,...
%  grainSurfaceSubscriptArray(:,3)-zBoundsNew(1)+1, 'b.'); 

line(xTopOfLoop*[1 1], [1 yTopOfLoop], zTopOfLoop*[1 1])

line(xBottomOfLoop*[1 1], [1 yBottomOfLoop], zBottomOfLoop*[1 1])

plot3(loopSubscriptArray(:,1), loopSubscriptArray(:,2),...
   loopSubscriptArray(:,3), 'r.'); 

plot3(curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3), 'b.')

plot3(xInterpolated, yToInterpolate, zToInterpolate, 'g.')



%% Cut up from centre line before cleaning.
for iSlice = 1:smallVolumeSize(3)
    % Find y positions of curve on slice.   
    tempIndexList = find(centreCurveVolume(:,:,iSlice));

    if ~isempty(tempIndexList)
        nIndex = length(tempIndexList); tempSubscriptArray = zeros(nIndex, 2);

        [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(smallVolumeSize(1:2), tempIndexList);
        
        [~, maxIndex] = max(tempSubscriptArray(:,2));
        
        % Get values in grain and edge volumes.
        volumeColumn = smallGrainVolume(tempSubscriptArray(maxIndex,1), (tempSubscriptArray(maxIndex,2)+1):end, iSlice);
        
        exteriorColumn = smallGrainExterior(tempSubscriptArray(maxIndex,1), (tempSubscriptArray(maxIndex,2)+1):end, iSlice);
        
        %Step up column and add into curve volume if air or edge
        for jPoint = 1:length(volumeColumn)
            
            if volumeColumn(jPoint) == 0 || exteriorColumn(jPoint)
                
                centreCurveVolume(tempSubscriptArray(maxIndex,1), tempSubscriptArray(maxIndex,2)+jPoint, iSlice) = 1;
                
                plot3(tempSubscriptArray(maxIndex,1), tempSubscriptArray(maxIndex,2)+jPoint, iSlice, 'k.')

            elseif volumeColumn(jPoint) == ENDOSPERM_INDEX || volumeColumn(jPoint) == GERM_INDEX
                %Stop once non edge endosperm or germ reached
               
               %plot3(tempSubscriptArray(maxIndex,1), tempSubscriptArray(maxIndex,2)+jPoint-1, iSlice, 'mx')
               
               break
            end
        end
    end
end


%% Check curve is 16 closed - on faces and edges
%%% Note diagonal holes in curve not closed - will allow [1 1; 1 0] on [0 1; 1 1]

%Close edges in 2D along each Y slice
for iSlice = 1:smallVolumeSize(3)
   
    % Check if exisiting crease points of crease on slice.
    tempIndexList = find(centreCurveVolume(:,:,iSlice));

    if ~isempty(tempIndexList)
        nIndex = length(tempIndexList); tempSubscriptArray = zeros(nIndex, 2);
        
        [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(smallVolumeSize(1:2), tempIndexList);
        
        if length(unique(tempSubscriptArray(:,2))) ~= max(tempSubscriptArray(:,2))
           % Should not occur afer interpolation...
            
           error('Missing Y values on slice') 
           %figure; plot(tempSubscriptArray(:,1), tempSubscriptArray(:,2),'.')
        end

        % Step up by rows.
        for jRow = 1:(max(tempSubscriptArray(:,2))-1)
        
            % Get subscripts on row
            xSubscriptList = find(centreCurveVolume(:,jRow,iSlice));

            % Test lateral conenctivity.
            if length(xSubscriptList) > 1
                if max(diff(xSubscriptList)) > 1 
                    error('Lateral connectivity broken.')
                end    
            end
                
            numberAbove = zeros(length(xSubscriptList),1);
            
            numberDiagonal = zeros(length(xSubscriptList),2);
            
            % Step along up subscripts checking if one has point above.
            for kStep = 1:length(xSubscriptList)
                
                if centreCurveVolume(xSubscriptList(kStep),jRow+1,iSlice)
                    
                    numberAbove(kStep) = 1;
                    
                    numberDiagonal(kStep,:) = NaN;
                else
                    %Check for diagonals
                    if centreCurveVolume(xSubscriptList(kStep)-1,jRow+1,iSlice)
                        numberDiagonal(kStep,1) = -1;
                    end
                    
                    if centreCurveVolume(xSubscriptList(kStep)+1,jRow+1,iSlice)
                        numberDiagonal(kStep,2) = 1;
                    end
                end
            end
            
            % If none with one above, need to fill diagonls
            if sum(numberAbove) == 0 
                
               if sum( abs(numberDiagonal(:))) 
                   
                   % If just one diagonal can do mathmatically.
                   if sum( abs(numberDiagonal(:))) == 1
                       
                       % Check which subscript diagonal belongs to.
                       for kStep = 1:length(xSubscriptList)
                           
                          if sum( abs(numberDiagonal(kStep,:)))
                              centreCurveVolume(xSubscriptList(kStep)+sum(numberDiagonal(kStep,:)), ...
                                  jRow,iSlice) = 1;
                              
                              plot3(xSubscriptList(kStep)+sum(numberDiagonal(kStep,:)), ...
                                  jRow, iSlice, 'k.')
                          end
                       end
                       
                   else
                      % This should only occur before a branch, 
                      % which would usually be flagged in next step. 
                       
                      error('Need to implement picking') 
                   end
                   
               else
                  % No diagonal connection to slice above, can occur on interpolated slices
                  
                  % Get subscripts on row above
                  xAboveList = find(centreCurveVolume(:,jRow+1,iSlice));

                  % Find closest subscripts on each row.
                  distance2Above = zeros(length(xSubscriptList),2);
                  
                  for kStep = 1:length(xSubscriptList)
                      
                      [distance2Above(kStep,1), distance2Above(kStep,2)] = ...
                          min( abs(xSubscriptList(kStep)-xAboveList ));
                  end
                  
                  [~, minIndex] = min(distance2Above(:,1));
                  
                  % Get points to fill.  
                  if xSubscriptList(minIndex) > xAboveList(distance2Above(minIndex,2))
                      
                      toFill = xAboveList(distance2Above(minIndex,2)):xSubscriptList(minIndex);
                  else
                      
                      toFill = xSubscriptList(minIndex):xAboveList(distance2Above(minIndex,2));
                  end
                  
                  if length(xAboveList) > 1 | length(xSubscriptList) > 1
                      pause(1);
                  end
                  
                  % Fill points on row
                  centreCurveVolume(toFill,jRow,iSlice) = 1;
                              
                  plot3(toFill, jRow, iSlice, 'r.')
               end
            end
        end
    end
end
    
%%% Test how is geodesic connectivity calculated - its 26
%temp = zeros(3,3,3,'logical');
%temp(1, :, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 1) = 1; temp(3, 3, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 2) = 1; temp(3, 3, 3) = 1;
%bwdistgeodesic(temp, 1,'quasi-euclidean')

%%% Add slight cut above end

%% Test geodesic embedding for unwrapping.
% Problem can occur if there are unlinked vertices, check by running fast marching once.
tempD = perform_fast_marching_mesh(aleuroneSurface.vertices, aleuroneSurface.faces, 50);

%%% If small number of unreachable faces, should be able to find all by testing
%%% one point randomly. Need to confirm.

% Remove points with inf distance.
if sum(isinf(tempD))
    vertexToRemove = find(isinf(tempD));

    facesToRemove = zeros(nFace,1);

    for iVertex = 1:length(vertexToRemove)
        vertexInFace = find(aleuroneSurface.faces == vertexToRemove(iVertex));

        [tempFaceIndex, ~] = ind2sub([nFace, 3], vertexInFace);

        facesToRemove(tempFaceIndex) = 1;
    end

    facesToRemove = find(facesToRemove);

    aleuroneSurface.faces(facesToRemove, :) = [];

    % Need to keep original vertex positions to relink faces
    initialVertexPositions = 1:nVertex;

    aleuroneSurface.vertices(vertexToRemove, :) = [];

    initialVertexPositions(vertexToRemove) = [];

    nVertex = size(aleuroneSurface.vertices,1);

    nFace = size(aleuroneSurface.faces,1);

    % Relink remaining faces.
    for iVertex = 1:nVertex

        aleuroneSurface.faces(aleuroneSurface.faces == initialVertexPositions(iVertex)) = iVertex;

    end
end

% Get distance map for all.
D = zeros(nVertex);

for iVertex = 1:nVertex
    
    D(:,iVertex) = perform_fast_marching_mesh(aleuroneSurface.vertices, aleuroneSurface.faces, iVertex);
    
end

sum(isinf(D(:)))

D = (D + D')/2;

J = eye(nVertex) - ones(nVertex)/nVertex;
W = -J*(D.^2)*J;

%%% Infs are sometimes produced in loop, need to resolve.

[U,S] = eig(W);
S = diag(S);
[S,I] = sort(S,'descend'); U = U(:,I);

% To plot.
vertexF = U(:,1:2)' .* repmat(sqrt(S(1:2)), [1 nVertex]);

icenter = 50;
irotate = 50;

vertexF = vertexF - repmat(vertexF(:,icenter), [1 nVertex]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];

figure;
plot_mesh(vertexF,aleuroneSurface.faces);
       
% To do: 
%        - Then apply deformation field to volumes
%        - Look down into volume to get thickness of aleurone
%        - Could also integrate intensity to see if this varies