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
%% Take mesh around entire grain then limit to aleurone exterior.

% Reduce patch size to save mesh construction time - keep as large as possible.
% If this value is below ~0.05, holes may appear in mesh if grain touches volume border.
% resizeRatio = 0.2;
% 
% % Do a strong morphological open to seperate along crease and remove floaters.
% smallGrainVolume = grainVolume;
% 
% % Volume should not be binarized before this point.
% smallGrainVolume = imresize3(uint8(smallGrainVolume), resizeRatio);
% 
% %%% Should test for presence of multiple objects/floaters here.
% 
% % Create isosurface on downsampled volume.
% fullGrainSurface = isosurface(smallGrainVolume, 0.5);
% 
% %voxelizedFullGrainSurface = inpolyhedron(fullGrainSurface, 1:voxelSize(1), ... 
% %   1:voxelSize(2), 1:voxelSize(3));
% 
% % Simplify patch to speed processing time later. 
% % This and resize ratio both need tuning to result in < 5000 vertices.
% patchReductionRatio = 0.02;
% 
% fullGrainSurface = reducepatch(fullGrainSurface,patchReductionRatio);
% 
% length(fullGrainSurface.vertices)
% 
% % Englarge surface back to original size.
% fullGrainSurface.vertices = round(fullGrainSurface.vertices/resizeRatio);
% 
% % Note that surface vertices are flipped relative to image subscripts.
% fullGrainSurface.vertices(:,[1 2]) = fullGrainSurface.vertices(:,[2 1]);
% 
% % Was testing mesh voxilization for debugging
% %voxelizedSmallGrainSurface = inpolyhedron(fullGrainSurface, 1:volumeSize(1), ... 
% %   1:volumeSize(2), 1:volumeSize(3));
% 
% % Truncate to lower bound of volume.
% fullGrainSurface.vertices(fullGrainSurface.vertices < 1) = 1;
%         
% % Truncate to upper bound.
% fullGrainSurface.vertices(fullGrainSurface.vertices(:,1) > volumeSize(1),1) = volumeSize(1);
% 
% fullGrainSurface.vertices(fullGrainSurface.vertices(:,2) > volumeSize(2),2) = volumeSize(2);
% 
% fullGrainSurface.vertices(fullGrainSurface.vertices(:,3) > volumeSize(3),3) = volumeSize(3);
% 
% nVertex = size(fullGrainSurface.vertices,1);
% 
% nFace = size(fullGrainSurface.faces,1);
% 
% % Get long axis of grain with PCA.
% grainAxisArray = pca(fullGrainSurface.vertices);
% 
% grainLongAxis = grainAxisArray(:,1)';
% 
% grainCreaseAxis = grainAxisArray(:,3)';
% 
% grainCenter = mean(fullGrainSurface.vertices);
% 
% % Plost test figures;
% figure; axis equal; hold on; axis off; set(gca, 'Clipping', 'off')
% 
% patch(fullGrainSurface);
% 
% plot3(fullGrainSurface.vertices(:,1), fullGrainSurface.vertices(:,2), fullGrainSurface.vertices(:,3), 'ro')
% 
% line(grainCenter(1)+[0 2000]*grainLongAxis(1), grainCenter(2)+[0 2000]*grainLongAxis(2), ...
%     grainCenter(3)+[0 2000]*grainLongAxis(3), 'color' ,'b')
% 
% line(grainCenter(1)+[0 1000]*grainAxisArray(1,2), grainCenter(2)+[0 1000]*grainAxisArray(2,2), ...
%     grainCenter(3)+[0 1000]*grainAxisArray(3,2), 'color' ,'g')
% 
% line(grainCenter(1)+[0 1000]*grainCreaseAxis(1), grainCenter(2)+[0 1000]*grainCreaseAxis(2), ...
%     grainCenter(3)+[0 1000]*grainCreaseAxis(3), 'color' ,'r')
% 
% title('Check there are no holes in the surface')



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
grainExterior = ~grainVolumeAligned;

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

clear grainExterior, clear grainIndexList, clear grainVolume

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

tempGrainVolume = grainVolumeAligned(xBoundsNew(1):xBoundsNew(2), ...
    yBoundsNew(1):yBoundsNew(2), zBoundsNew(1):zBoundsNew(2));

tempSize = size(tempGrainVolume);



% Find loop start on top of grain. Firstly step in from top of grain
zTopOfLoop = zRange(1)-zBoundsNew(1)+1 + 10;

% Then take centre X position.
sliceIndexList = find(tempGrainVolume(:,:,zTopOfLoop));

nIndex = length(sliceIndexList); sliceSubscriptArray = zeros(nIndex, 3);

[sliceSubscriptArray(:,1), sliceSubscriptArray(:,2), sliceSubscriptArray(:,3)] = ...
    ind2sub(tempSize, sliceIndexList);

xTopOfLoop = round( mean(sliceSubscriptArray(:,1)));

% Then take lowest Y in line.
tempIndex = find(sliceSubscriptArray(:,1) == xTopOfLoop);

yTopOfLoop = min(sliceSubscriptArray(tempIndex,2))-1;

% Repeat for bottom of loop.
zBottomOfLoop = zRange(2)-zBoundsNew(1)+1 - 10;

sliceIndexList = find(tempGrainVolume(:,:,zBottomOfLoop));

nIndex = length(sliceIndexList); sliceSubscriptArray = zeros(nIndex, 3);

[sliceSubscriptArray(:,1), sliceSubscriptArray(:,2), sliceSubscriptArray(:,3)] = ...
    ind2sub(tempSize, sliceIndexList);

xBottomOfLoop = round( mean(sliceSubscriptArray(:,1)));

tempIndex = find(sliceSubscriptArray(:,1) == xBottomOfLoop);

yBottomOfLoop = min(sliceSubscriptArray(tempIndex,2))-1;



% Calc distance map from plane underneath grain.
%basePlaneVolume = zeros(tempSize, 'logical');

%basePlaneVolume(:,1,:) = 1;

% Make mask image for grain.
grainMask = logical(tempGrainVolume);

% Step through and cut above centre line;
for iSlice = 1:tempSize(3)
    
    if ~isnan( centreLineHeightProfile(iSlice+zBoundsNew(1)-1))
        
        cutAbove = centreLineHeightProfile(iSlice+zBoundsNew(1)-1)-yBoundsNew(1)+1+50;
        
        if cutAbove < tempSize(2)
            
            grainMask(:, cutAbove:end, iSlice) = 1;
        end 
    else
        % Cut above centre of mass of pixels.
        sliceIndexList = find(tempGrainVolume(:,:,iSlice));

        [~, ySubscriptList, ] = ind2sub(tempSize, sliceIndexList); 
        
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

% Graydist takes a long time.
tic
distanceFromTop = graydist(distanceFromGrain, sub2ind(tempSize, ...
    xTopOfLoop, yTopOfLoop, zTopOfLoop), 'quasi-euclidean');
toc

tic
distanceFromBottom = graydist(distanceFromGrain, sub2ind(tempSize, ...
    xBottomOfLoop, yBottomOfLoop, zBottomOfLoop),'quasi-euclidean');
toc

dMap = distanceFromTop + distanceFromBottom;

clear distanceFromTop, clear distanceFromBottom

% Round to lower precision to prevent floating point errors.
%%% Note, changed from 8 to 32 as loop doesn't reach end.
dMap = round(dMap * 64)/64;

dMap(isnan(dMap)) = Inf;



% Reginal minima should define connection.
loopVolume = imregionalmin(dMap);

clear dMap



% Take largest connected region of loop
tempCC = bwconncomp(loopVolume, 26);

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
    loopVolume(tempStats(iRegion).PixelIdxList) = 0;
end

loopIndexList = find(loopVolume);

nIndex = length(loopIndexList); loopSubscriptArray = zeros(nIndex, 3);

[loopSubscriptArray(:,1), loopSubscriptArray(:,2), loopSubscriptArray(:,3)] = ...
    ind2sub(tempSize, loopIndexList);

zTopOfLoop = min(loopSubscriptArray(:,3));

zBottomOfLoop = max(loopSubscriptArray(:,3));



% Test plot.
figure; hold on; axis equal; set(gca, 'Clipping', 'off')

line(xTopOfLoop*[1 1], [1 yTopOfLoop], zTopOfLoop*[1 1])

line(xBottomOfLoop*[1 1], [1 yBottomOfLoop], zBottomOfLoop*[1 1])

plot3(loopSubscriptArray(:,1), loopSubscriptArray(:,2),...
   loopSubscriptArray(:,3), 'r.'); 

%%% Expand on this, need to create connected line between slices
%%% Step along, taking highest connected points on each slice
% Take highest point above centre of loop at each slice.
% for iSlice = topOfLoopZ:bottomOfLoopZ
%     
%     loopSlice = loopVolume(:,:,iSlice);
%     
%     % Check if more than one coordinate.
%     if sum(loopSlice(:)) > 1
%         % Take average of coordinates.
%         tempIndexList = find(loopSlice);
% 
%         nIndex = length(tempIndexList); tempSubscriptArray = zeros(nIndex, 2);
% 
%         [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(tempSize(1:2), tempIndexList);
%         
%         averagePosition = mean(tempSubscriptArray);
%         
%         % Find closest point on exisiting loop to average position
%         [~, minIndex] = min( sqrt((averagePosition(1) - tempSubscriptArray(:,1)).^2 + ... 
%                 (averagePosition(2) - tempSubscriptArray(:,2)).^2));
%                 
%         averagePosition = tempSubscriptArray(minIndex,:);
%         
%         % Raise average position up until it is just below volume
% %         while grainMask(averagePosition(1),averagePosition(2)+1,iSlice) == 0
% %             averagePosition(2) = averagePosition(2)+1;
% %             
% %             if averagePosition(2) == tempSize(2)
% %                error('Midline has risen to top of volume') 
% %             end
% %         end
% 
%         % Clear loop slice and place single index.
%         loopSlice(tempIndexList) = 0;
%         
%         loopSlice(averagePosition(1), averagePosition(2)) = 1;
%         
%         loopVolume(:,:,iSlice) = loopSlice;
%     end
% end
%
%
%
% loopIndexList = find(loopVolume);
% 
% nIndex = length(loopIndexList); loopSubscriptArray = zeros(nIndex, 3);
% 
% [loopSubscriptArray(:,1), loopSubscriptArray(:,2), loopSubscriptArray(:,3)] = ...
%     ind2sub(tempSize, loopIndexList);
% 
% plot3(loopSubscriptArray(:,1), loopSubscriptArray(:,2),...
%    loopSubscriptArray(:,3), 'bx'); 

%% Now find centre-curve by slice.

centreCurveVolume = zeros(tempSize, 'logical');

% Transform distance from grain so centrelines are emphasized.
distanceNearGrain = distanceFromGrain;

distanceNearGrain(grainMask) = 0;

distanceNearGrain = -log(distanceNearGrain);

distanceNearGrain = distanceNearGrain - min(distanceNearGrain(:));

distanceNearGrain(grainMask) = Inf;

basePlane = zeros(tempSize(1:2), 'logical');

basePlane(:,1) = 1;

missingSlices = zeros(tempSize(3),1);

for iSlice = 1:tempSize(3)
   
    % Check if exisiting crease points of crease on slice.
    indexList = find(loopVolume(:,:,iSlice));
    
    if ~isempty(indexList)
        
        % Calculate distance maps to find minimum path from top of crease to base plane.
        distanceFromCrease = graydist(distanceNearGrain(:,:,iSlice), loopVolume(:,:,iSlice), 'quasi-euclidean');
        
        distanceFromBase = graydist(distanceNearGrain(:,:,iSlice), basePlane, 'quasi-euclidean');
        
        dMap = distanceFromCrease + distanceFromBase;
        
        dMap = round(dMap * 64)/64;

        dMap(isnan(dMap)) = Inf;
        
        % Regional minima is shortest path
        curveSlice = imregionalmin(dMap);

        % If regions are not linked, everything will be Inf
        curveSlice(curveSlice & isinf(dMap)) = 0;

        %temp = curveSlice + grainMask(:,:,iSlice);
        %figure; imshow(temp)
        

        
        % Test that line links top and bottom.
        planeInCurve = sum(curveSlice(:) & basePlane(:));
        
        loopInCurve = sum(sum(curveSlice & loopVolume(:,:,iSlice)));

        if planeInCurve & loopInCurve
            centreCurveVolume(:,:,iSlice) = curveSlice;

            curveIndexList = find(curveSlice);

            nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

            [curveSubscriptArray(:,1), curveSubscriptArray(:,2)] = ...
               ind2sub(tempSize(1:2), curveIndexList);

           plot3(curveSubscriptArray(:,1), curveSubscriptArray(:,2),...
                iSlice*ones(nIndex,1), 'b.'); 
        else
           missingSlices(iSlice) = 1; 
        end
    end
end

%% Interpolate missing slices.
missingSliceIndexList = find(missingSlices);

sliceSteps = find(diff(missingSliceIndexList) > 1);

nMissingSlices = length(missingSliceIndexList);

% Get points to interpolate.
yToInterpolate = zeros(nMissingSlices, tempSize(2))*NaN;

zToInterpolate = zeros(nMissingSlices, tempSize(2))*NaN;

for iSlice = 1:nMissingSlices
    % Find y position of loop on slice   
    tempIndexList = find(loopVolume(:,:,missingSliceIndexList(iSlice)));

    nIndex = length(tempIndexList); tempSubscriptArray = zeros(nIndex, 2);

    [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(tempSize(1:2), tempIndexList);
    
    % Only allocate up to highest Y on loop.
    yToInterpolate(iSlice,1:max(tempSubscriptArray(:,2))) = 1:max(tempSubscriptArray(:,2));
    
    zToInterpolate(iSlice,:) = missingSliceIndexList(iSlice);
end

% Remove excess interpolation points.
zToInterpolate(isnan(yToInterpolate)) = [];

yToInterpolate(isnan(yToInterpolate)) = [];



curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(tempSize, curveIndexList);

xInterpolated = griddata(curveSubscriptArray(:,2), curveSubscriptArray(:,3), curveSubscriptArray(:,1), ...
    yToInterpolate, zToInterpolate, 'linear');

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

%%% Check it is closed
closureTest = convn(centreCurveVolume, STREL_18_CONNECTED.Neighborhood);

figure; hist(closureTest(closureTest ~= 0));

%Test each point for 16 closure, if not, draw line to adjacent point on x
%slice

%%% Add slight cut in and above end


%% Split alonge crease by slice. - probably wont use.
for iSlice = 1900; zRangeAleurone(1):20:zRangeAleurone(2)
    if all( ~isnan( boundsLineBySlice(iSlice,:)))
        % Set outside of bounds to one on slice.
        exteriorVolume(1:boundsLineBySlice(iSlice,1),:,iSlice) = 1;
        
        exteriorVolume(boundsLineBySlice(iSlice,2):end,:,iSlice) = 1;
%         
%         maskIm = ~exteriorVolume(:,:,iSlice);
%         DMask = bwdist(exteriorVolume(:,:,iSlice), 'quasi-euclidean');
%         
%         % Calculate line from base of volume to top of centre line.
%         bwBase = zeros(volumeSize(1:2)); bwBase(centreLine(iSlice),1) = 1;
%         D1 = bwdist(bwBase, 'quasi-euclidean');
%         
%         topInd = find(exteriorVolume(centreLine(iSlice),:,iSlice) == 0);
%         bwTop = zeros(volumeSize(1:2)); bwTop(centreLine(iSlice),topInd(end)) = 1;
%         D2 = bwdist(bwTop, 'quasi-euclidean');
%         DMap = D1 + D2 + DMask;
%         
%         paths = imregionalmin(DMap);
%         P = imoverlay(maskIm, paths, [.5 .5 .5]);
% 
%         figure;
%         imshow(P);

        % Calculate distance map to use for weight later on.
        %mask = ~exteriorVolume(:,:,iSlice);
        %start = 

        tempImage = bwdistgeodesic(~exteriorVolume(:,:,iSlice), 'quasi-euclidean');
        
        maxHeightOnSlice = max(creaseProfileBySlice(iSlice,:,2));
        
        sliceCentre = zeros(maxHeightOnSlice,1)*NaN;
   
        % Range of X to include in center calculation.
        xIncluded = boundsLineBySlice(iSlice,1):boundsLineBySlice(iSlice,2); 
        
        figure; imshow(tempImage); hold on;

        % Step up y in rows of x.
        for jRow = 1:maxHeightOnSlice
            % Remove X if above height for column.    
            temp = find(creaseProfileBySlice(iSlice,xIncluded,2) < jRow);
            
            xIncluded(temp) = [];

            % Find x values corresponding to air.
            temp = find(grainVolumeAligned(xIncluded,jRow,iSlice) == 0);
            
            if length(temp) > 3
                % Take average centreline given distance weighting
                sliceCentre(jRow) = wmean(xIncluded(temp), tempImage(xIncluded(temp),jRow)');
            end
            
            %plot(jRow, xIncluded, 'g.')
        end
        
        %Interp any missing values.
        tempMissing = find( isnan(sliceCentre));

        tempHave = find( ~isnan(sliceCentre));

        sliceCentre(tempMissing) = interp1(tempHave, sliceCentre(tempHave), ...
            tempMissing, 'linear');

        % Smooth.
        tempHave = find( ~isnan(sliceCentre));

        sliceCentre(tempHave) = round( smooth(sliceCentre(tempHave)));
        
        sliceCentre(isnan(sliceCentre)) = [];
        
        sliceCentre(end+1:end+10) = sliceCentre(end);
        
        plot(1:length(sliceCentre), sliceCentre, 'r.')
     end
end

%% Snap each surface vertices onto the surface voxels and test if aleurone.
aleuroneSurface = fullGrainSurface;

vertexToRemove = zeros(nVertex,1);

facesToRemove = zeros(nFace,1);

voxelRange = 25;

for iVertex = 1:nVertex
    surfaceVertex = aleuroneSurface.vertices(iVertex,:);
    
    % If vertex is not already a surface voxel.
    if ~grainExterior(surfaceVertex(1), surfaceVertex(2), surfaceVertex(3))
        % Find closest voxel and update vertex with this.

%         if iVertex == 59
%         % Calculation is slow due to number of voxels to test
%             [~, closestVoxelIndex] = min( sqrt( (grainSurfaceSubscriptArray(:,1)-surfaceVertex(1)).^2 + ...
%                 (grainSurfaceSubscriptArray(:,2)-surfaceVertex(2)).^2 + ...
%                 (grainSurfaceSubscriptArray(:,3)-surfaceVertex(3)).^2 ) );
%             grainSurfaceSubscriptArray(closestVoxelIndex,:)
%         end
        
        %%% ADD FUNCTION for volume test with and add above
        % Check subvolume range is inside main volume, truncate if not.
        subVolumeRange = [surfaceVertex-voxelRange; surfaceVertex+voxelRange];
        
        % Keep track of voxel offset from zero for later.
        offset = (voxelRange + 1)*[1 1 1];
        
        % Change offset if lower volume bound is below 1.
        if subVolumeRange(1,1) < 1; offset(1) = offset(1) + subVolumeRange(1,1) - 1; end
        
        if subVolumeRange(1,2) < 1; offset(2) = offset(2) + subVolumeRange(1,2) - 1; end
        
        if subVolumeRange(1,3) < 1; offset(3) = offset(3) + subVolumeRange(1,3) - 1; end
        
        % Truncate to lower bound.
        subVolumeRange(subVolumeRange < 1) = 1;
        
        % Truncate to upper bound.
        if subVolumeRange(2,1) > volumeSize(1); subVolumeRange(2,1) = volumeSize(1); end
        
        if subVolumeRange(2,2) > volumeSize(2); subVolumeRange(2,2) =  volumeSize(2); end
        
        if subVolumeRange(2,3) > volumeSize(3); subVolumeRange(2,3) = volumeSize(3); end
        
        subVolumeDimensions = subVolumeRange(2,:) - subVolumeRange(1,:) + 1;
        
        % Get subvolume near vertex to speed up.
        tempSubvolume = grainExterior(subVolumeRange(1,1):subVolumeRange(2,1), ...
            subVolumeRange(1,2):subVolumeRange(2,2), subVolumeRange(1,3):subVolumeRange(2,3));
        
        % Get indices of surface.
        tempSurfaceIndexList = find(tempSubvolume);

        % If some surface indices in subvolume, calculate using this.
        if ~isempty(tempSurfaceIndexList)

            tempNo = length(tempSurfaceIndexList);

            % Convert indices to subvolumes.
            tempSubscriptArray = zeros(tempNo, 3);
            
            [tempSubscriptArray(:,1), tempSubscriptArray(:,2), tempSubscriptArray(:,3)] = ...
                ind2sub(subVolumeDimensions, tempSurfaceIndexList);
            
            % Offset back to original coordinates.
            tempSubscriptArray = [tempSubscriptArray(:,1) + surfaceVertex(1) - offset(1), ...
                tempSubscriptArray(:,2) + surfaceVertex(2) - offset(2), ...
                tempSubscriptArray(:,3) + surfaceVertex(3) - offset(3)];
            
            % Find closest voxel and update vertex with this.
            [~, closestVoxelIndex] = min( sqrt( (tempSubscriptArray(:,1)-surfaceVertex(1)).^2 + ...
                (tempSubscriptArray(:,2)-surfaceVertex(2)).^2 + ...
                (tempSubscriptArray(:,3)-surfaceVertex(3)).^2 ) );
            
            surfaceVertex = tempSubscriptArray(closestVoxelIndex,:);

            aleuroneSurface.vertices(iVertex,:) = surfaceVertex;
            
        else
            
           error('No surface voxels in subvolume') 
           
        end
    end
    
    
    
    % Test if vertex is part of aleurone and mark to remove if not.
    voxelMaterialIndex = grainVolume(surfaceVertex(1), surfaceVertex(2), surfaceVertex(3));
    
    if ALEURONE_INDEX ~= voxelMaterialIndex
        vertexToRemove(iVertex) = 1;
        
        % Also mark faces to remove.
        vertexInFace = find(aleuroneSurface.faces == iVertex);
        
        [tempFaceIndex, ~] = ind2sub([nFace, 3], vertexInFace);
        
        facesToRemove(tempFaceIndex) = 1;
    end
end

%%% ADD FUNCTION for vertex/face removal and use below
% Remove faces, vertices, and update numbers.
facesToRemove = find(facesToRemove);

aleuroneSurface.faces(facesToRemove, :) = [];

vertexToRemove = find(vertexToRemove);

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

aleuroneAxisArray = pca(aleuroneSurface.vertices);

aleuroneCreaseAxis = aleuroneAxisArray(:,3)';



%% Plot test figures.
figure; axis equal; hold on; set(gca, 'Clipping', 'off')

iToPlot = find(aleuroneSurfaceRotated.vertices(:,2) < 75);

%patch(aleuroneSurfaceRotated); 

plot3(aleuroneSurfaceRotated.vertices(iToPlot,1), aleuroneSurfaceRotated.vertices(iToPlot,2), ...
   aleuroneSurfaceRotated.vertices(iToPlot,3), 'r.');

%split into two clusters
idx = kmeans(aleuroneSurfaceRotated.vertices(iToPlot,[1 2]),2);

cluster1index = find(idx == 1);

cluster2index = find(idx == 2);

plot3(aleuroneSurfaceRotated.vertices(iToPlot(cluster1index),1), aleuroneSurfaceRotated.vertices(iToPlot(cluster1index),2), ...
    aleuroneSurfaceRotated.vertices(iToPlot(cluster1index),3), 'r.');

plot3(aleuroneSurfaceRotated.vertices(iToPlot(cluster2index),1), aleuroneSurfaceRotated.vertices(iToPlot(cluster2index),2), ...
    aleuroneSurfaceRotated.vertices(iToPlot(cluster2index),3), 'b.');

% Histogram provides some seperation.
figure; hist(aleuroneSurfaceRotated.vertices(iToPlot,1),200)

figure; plot(aleuroneSurfaceRotated.vertices(iToPlot,1), aleuroneSurfaceRotated.vertices(iToPlot,2), '.')



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