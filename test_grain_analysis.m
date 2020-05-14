% To do
    %Normalize intensity between multiple grains based on germ histogram

% Test load and plot of grain stack.
clc; clear; close all

% Other: Om_1_7_test
% labelDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Grain LU/Data/OB6_test/Labels';
% greyDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Grain LU/Data/OB6_test/Images';

labelDirectory = 'C:\Users\Admin\Lund University\Nick Sirijovski - Gavin\Data_for_ML_Test\OB6_test\Labels';
greyDirectory = 'C:\Users\Admin\Lund University\Nick Sirijovski - Gavin\Data_for_ML_Test\OB6_test\Images';

% Loading tiff stack takes a while.
grainVolume = loadtiffstack(labelDirectory, 1);

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

VOXEL_SIZE = 4;

% Flags
calculateArea = 0; % Very time consuming.

testPlotNormals = 0;

maxAleuroneThickness = 75; %In voxels - previous limit at 20

blockThickness = 20/VOXEL_SIZE;

depthToCalculate = 200/VOXEL_SIZE;
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

grainVolumeAligned = uint8( affine_transform_full(single(grainVolume), M1*M2*M3*M4, 5));

%%% Note small holes appear after affine transform, close fixes these. 
%%% Probably bug with nearest neighbour interp I added to affine transform c file
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

% Extrapolate polygon to top and bottom of volume


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
%%% Note, changed from 8 to 64 as loop doesn't reach end,
%%% Causes some floaters, but no gaps. Not guaranteed to work for all...
dMap = round(dMap * 64)/64;

dMap(isnan(dMap)) = Inf;

% Reginal minima should define connection.
loopVolume = imregionalmin(dMap);

clear dMap

%% Take top of connected curve.
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

clear distanceFromEnd

%Clear loop and add replace with thin values.
loopVolume(:) = 0;

tempIndexList = sub2ind(smallVolumeSize, coordinatesOfLoop(:,1), coordinatesOfLoop(:,2), coordinatesOfLoop(:,3));

loopVolume(tempIndexList) = 1;

loopIndexList = find(loopVolume);

nIndex = length(loopIndexList); loopSubscriptArray = zeros(nIndex, 3);

[loopSubscriptArray(:,1), loopSubscriptArray(:,2), loopSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, loopIndexList);

%% Now find centre-curve by slice.

% Transform distance from grain so centrelines are emphasized.
% Removed, seems to work better with linear distance
%distanceFromGrain = distanceFromGrain.^2;

distanceFromGrain(grainMask) = 0;

centreCurveVolume = zeros(smallVolumeSize, 'logical');

startIndex = sub2ind(smallVolumeSize(1:2), roughXCentre, 1);

for iSlice = [zTopOfLoop:5:(zBottomOfLoop-5+1) zBottomOfLoop]
   
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
                    (curveLine(jStep+1,2)-curveLine(jStep,2))^2) < 10

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

        if sum(tempMask(tempIndexList))
%                 figure; imshow(grainMask(:,:,iSlice));
% %%%              figure; imshow(distanceFromBase/50000);  
%                 hold on; plot(curveLine(:,2), curveLine(:,1),'r')
%                 plot(tempSubscriptArray(2), tempSubscriptArray(1),'gx')
%                  plot(curveLineFull(:,2), curveLineFull(:,1),'m.')
        end

        % Removed test for multiple connections to base plane, fast marching prevents.

        % Test that line links top and bottom.
        if curveSlice(roughXCentre,1) && ...
                curveSlice(tempSubscriptArray(1), tempSubscriptArray(2))
            % Removed thining and multi-part removal code, fast marching prevents.

            % Add slice of curve and loop to curve volume.
            centreCurveVolume(:,:,iSlice) = curveSlice | loopVolume(:,:,iSlice);

            % Vertical continuity test removed, fast marching should fix.
        else
           % Start and end don't connect, just add exisiting loop. 
           %centreCurveVolume(:,:,iSlice) = loopVolume(:,:,iSlice);

        end
    end
end

clear distanceFromGrain, clear loopVolume

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
toInterp = 1:10:size(xArray,2);

% Add final point if not included
if toInterp(end) ~= size(xArray,2)
    toInterp = [toInterp size(xArray,2)];
end

tic
curveNrb = nrbloft_nan(xArray(:, toInterp), yArray(:, toInterp), zArray(:, toInterp), 2);
toc

%% Put interpolated values into curve volume.
%zToInterp = zBottomOfLoop-zTopOfLoop+1;
zToInterp = max(max(zArray(:, toInterp)))-zTopOfLoop+1;

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

% Fill up to max X on each point.
centreCurveVolume(:) = 0;

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
centreCurveVolume = logical(centreCurveVolume - imerode(centreCurveVolume, STREL_26_CONNECTED));

% Remove ends from volume.
centreCurveVolume(:,:,1:(zTopOfLoop-1)) = 0;

centreCurveVolume(:,:,(zToInterp+1):smallVolumeSize(3)) = 0;

% Remove points above curve from volume
for iSlice = zTopOfLoop:zToInterp
    
   for jRow = (maxYValues(iSlice)+1):smallVolumeSize(2)

      centreCurveVolume(:, jRow, iSlice) = 0;

   end
end

%% Cut up from centre line and ends 
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

for iSlice = 1:smallVolumeSize(3)
    % Find y positions of curve on slice.   
    tempIndexList = find(centreCurveVolume(:,:,iSlice));

    if ~isempty(tempIndexList)
        nIndex = length(tempIndexList); tempSubscriptArray = zeros(nIndex, 2);

        [tempSubscriptArray(:,1), tempSubscriptArray(:,2)] = ind2sub(smallVolumeSize(1:2), tempIndexList);
        
        [~, maxIndexList] = max(tempSubscriptArray(:,2));
        
        maxIndexList = find(tempSubscriptArray(:,2) == tempSubscriptArray(maxIndexList,2));
        
        % If multiple max y values on one slice cut up from each to ensure connectivity.
        
        for jIndex = 1:length(maxIndexList)
            % Get values in grain and edge volumes.
            volumeColumn = smallGrainVolume(tempSubscriptArray(maxIndexList(jIndex),1), ...
                (tempSubscriptArray(maxIndexList(jIndex),2)+1):end, iSlice);

            exteriorColumn = smallGrainExterior(tempSubscriptArray(maxIndexList(jIndex),1), ...
                (tempSubscriptArray(maxIndexList(jIndex),2)+1):end, iSlice);

            %Step up column and add into curve volume if air or edge
            for jPoint = 1:length(volumeColumn)

                if volumeColumn(jPoint) == 0 || exteriorColumn(jPoint) 
                    %|| volumeColumn(jPoint) == ALEURONE_INDEX

                    centreCurveVolume(tempSubscriptArray(maxIndexList(jIndex),1), ...
                        tempSubscriptArray(maxIndexList(jIndex),2)+jPoint, iSlice) = 1;

%                     plot3(tempSubscriptArray(maxIndexList(jIndex),1),...
%                         tempSubscriptArray(maxIndexList(jIndex),2)+jPoint, iSlice, 'k.')

                % Stop once non- edge or air reach to prevent disconnection.
                % Only do if in main portion of grain, borders should extend up
                elseif iSlice > zTopOfLoop && iSlice < zToInterp
                    %if volumeColumn(jPoint) == ENDOSPERM_INDEX || volumeColumn(jPoint) == GERM_INDEX
                    
                    % plot3(tempSubscriptArray(maxIndexList(jIndex),1), tempSubscriptArray(maxIndexList(jIndex),2)+jPoint-1, iSlice, 'mx')

                   break
                end
            end
        end
    end
end

clear smallGrainExterior, clear smallGrainVolume
%% Test plot.
figure; hold on; axis equal; set(gca, 'Clipping', 'off')

line(xTopOfLoop*[1 1], [1 yTopOfLoop], zTopOfLoop*[1 1])

line(xBottomOfLoop*[1 1], [1 yBottomOfLoop], zBottomOfLoop*[1 1])

plot3(loopSubscriptArray(:,1), loopSubscriptArray(:,2),...
   loopSubscriptArray(:,3), 'r.'); 

plot3(curveSubscriptArray(:,1), curveSubscriptArray(:,2), ...
    curveSubscriptArray(:,3), 'b.')

% plot3(curveSubscriptArray(1:10:end,1), curveSubscriptArray(1:10:end,2), ...
%     curveSubscriptArray(1:10:end,3), 'b.')

% 1st dimension is Y, 2nd dimension is Z
%nrbplot(curveNrb, [10, 50]);

%plot3(curveX, curveY, curveZ, 'm.')

%% Place curve cut back into full volume
curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, curveIndexList);

% Correct subscripts back to original coordinates
curveSubscriptArray(:,1) = curveSubscriptArray(:,1) + xBoundsNew(1) - 1;

curveSubscriptArray(:,2) = curveSubscriptArray(:,2) + yBoundsNew(1) - 1;

curveSubscriptArray(:,3) = curveSubscriptArray(:,3) + zBoundsNew(1) - 1;

curveIndexList = sub2ind(volumeSize, curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3));

grainCutVolume = zeros(volumeSize, 'logical');

grainCutVolume(tempIndexList(grainExterior(tempIndexList) == 1)) = 1;

% Code curve into grain exterior as 2 to distinguish it later.
%grainExterior = uint8(grainExterior);
grainExterior(tempIndexList(grainExterior(tempIndexList) == 1)) = 0;

figure; imshow(sum(grainExterior(:,:,1623),3)/2)

clear centreCurveVolume 

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
germExterior = (grainExterior | grainCutVolume) & (grainVolumeAligned == GERM_INDEX);

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
%% Get exterior surface except for crease
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

for iRegion = 1:2
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
%% Check that aleurone and endosperm form continous region
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
combinedEdge = imdilate((grainExterior | grainCutVolume) & ~combinedExterior, STREL_18_CONNECTED);

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

clear combinedExterior, clear combinedEdge, clear grainCutVolume

%% Calculate surface area of aleurone exterior and interior
%https://se.mathworks.com/matlabcentral/answers/93023-is-there-a-matlab-function-that-can-compute-the-area-of-my-patch
% Matched Amira results closely for bee eye, should test for grain

if calculateArea
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
end
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

% for 300, 500 gives 12 9 points 
% for 20, 50, gives 343 1314 points 
% for 50, 100, gives 121 309 points 
% for 50, 75, gives 121 565 points 
% for 30, 75, gives 216 562 points  
% for 5, 20, gives 2072 8742 points  

% Run time:
% 60s per point on Mac
% 17s on Windows, 
    % 50-100s in parfor w/ 6 workers (picks up speed to 20-30 at end) avg 12s
    % 25-40s in parfor w/ 4 workers avg 8s

% For test, seems good to get dense on border and moderate on surface

edgeDistance = 20; %10

surfaceDistance = 50; %50

%Not, sparse distance must be less then others for code to work
sparsePointsDistance = 3; %3

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

%             testSparse = sum( sqrt( (aleuroneSurfaceSubscriptArray(sparsePointsChoosen,1) - pointChoosen(1)).^2 + ...
%                 (aleuroneSurfaceSubscriptArray(sparsePointsChoosen,2) - pointChoosen(2)).^2 + ...
%                 (aleuroneSurfaceSubscriptArray(sparsePointsChoosen,3) - pointChoosen(3)).^2) < sparsePointsDistance) == 0;
%             if testSparse
                 sparsePointsChoosen(sparsePointsToChoose(pointsOnLayer(ind))) = 1;
%             end

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

greyVolumeAligned = uint8( affine_transform_full(single(greyVolume), M1*M2*M3*M4, 1));

% Affine transform with linear interp. does not seem to have small holes problem.

clear greyVolume
%% Calulate distance between map points

%%% Test how is geodesic connectivity calculated - its 26
%temp = zeros(3,3,3,'logical');
%temp(1, :, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 1) = 1; temp(3, 3, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 2) = 1; temp(3, 3, 3) = 1;
%bwdistgeodesic(temp, 1,'quasi-euclidean')

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

% Set up for normal calculation
normalRadius = 50;

grid3D.nx = volumeSize(1); grid3D.ny = volumeSize(2); grid3D.nz = volumeSize(3);
grid3D.minBound = [1 1 1]';
grid3D.maxBound = volumeSize';

% Cut label volume so normals see seperated sides
    % Can cause problems with normal picking or step back if not
    % Note that no test points should be on cut as they come from grain exterior
grainVolumeAlignedCut = grainVolumeAligned;

grainVolumeAlignedCut(curveIndexList) = 0;

% Calculate normals to get thickness
%parfor iPoint = 1:nPoints %
for iPoint = 1154; 1261:nPoints
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
        if isempty(forwardAirIndex) & ~isempty(backwardAirIndex)
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

                
%                 figure; imshow(grainVolumeAlignedCut(:,:,currentSubscript(3))*100)
%                 hold on
%                 
%                 line( [0 tempNormal(2)*200]+currentSubscript(2),  ...
%                     [0 tempNormal(1)*200]+currentSubscript(1),...
%                         'color', 'r')
%                     
%                 plot(currentSubscript(2), currentSubscript(1) , 'c*')
% 
%                 plot(combinedSurfaceSubscripts(indexListInRange,2),...
%                     combinedSurfaceSubscripts(indexListInRange,1), '.')
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

            if indexList(1) == tempIndexList(1)
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
            else
               error('%i %i - Start indexes do not match', iPoint, jSubscript) 
            end
            
            % Find intersects
            % Using interior intersects some extra skimming problems
            lineIDs = grainVolumeAlignedCut(indexList);
            
            interiorList = aleuroneInterior(indexList);
            
            interiorPoints = find(interiorList);

            %[lineIDs(1:50) interiorList(1:50)]
            
            % Test there is an intersect
            if isempty(interiorPoints)
                if any(lineIDs == ENDOSPERM_INDEX) 
                   % This should not really occur...
                   
                   % For now, just add an interior intersect on al voxel before first endosperm
                   endospermInds = find(lineIDs == ENDOSPERM_INDEX);
                   
                   endospermInds = endospermInds(1);
                   
                   % Check voxel before is aleurone
                   if lineIDs(endospermInds - 1) == ALEURONE_INDEX
                       % If so, add in as (only interior intersect)
                       interiorPoints = endospermInds - 1;
                       
                       interiorList(interiorPoints) = 1;
                       
                       warning('%i %i - Added missing interior', iPoint, jSubscript)
                   else
                       error('%i %i - Hit endosperm without interior intersect', iPoint, jSubscript)
                   end
               else
                   % This should never occur 
                   error('%i %i - Not interior intersect or endosperm hit', iPoint, jSubscript) 
                end
            end
            
            intersectsCleared = 0;
            
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
                        error('%i %i - Unexpected materials before intersect', iPoint, jSubscript)
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

save(sprintf('C:\\Users\\Admin\\Documents\\MATLAB\\Temp_data\\%s_%i_%i_%i_%i_wEndo_distOnFull_voxelFractions_updateSkimming_splitCalcs', 'distanceMatrix', ...
        edgeDistance, surfaceDistance, sparsePointsDistance, normalRadius), ...
    'edgeDistance', 'surfaceDistance', 'sparsePointsDistance', 'normalRadius',...    
    'distanceMatrix', 'subscriptsToInterpolate', 'interpolatedIdentity',... 
    'distanceMatrixSparse', 'subscriptsForSparse', 'sparseIdentity',...
    'normalByPoint', 'internalIntersectByPoint', 'thicknessByPoint', 'averageIntensityByPoint',...
    'normalForSparse', 'internalIntersectForSparse', 'thicknessForSparse', 'averageIntensityForSparse',...
    'pointToSparseLinks', 'indexsToUseByPoint', '-v7.3');

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
    
%% Calaculate integrals under surface.

% Try to make these all integer values given Voxel thickenss
numberOfBlocks = depthToCalculate/blockThickness;

pointIntensityProfile = zeros(nPoints, numberOfBlocks)*NaN;

sparseIntensityProfile = zeros(nSparsePoints, numberOfBlocks)*NaN;

% Some catches included to correct for some points being lost if normal
% skims border of other structure. - but some have to be lost...

% Step through main points
for iPoint = 1:(nPoints + nSparsePoints)
    
    if iPoint <= nPoints
        tempNormal = normalByPoint(iPoint,:);
        
        pointIdentity = interpolatedIdentity(iPoint);
    else
        tempNormal = normalForSparse(iPoint-nPoints,:);
        
        pointIdentity = sparseIdentity(iPoint-nPoints);
    end
    
    if ~any(isnan(tempNormal))
        
        if iPoint <= nPoints
            %Check if aleurone or endosperm
            if pointIdentity
                %Aleurone, extend line in from interior surface
                startPoint = internalIntersectByPoint(iPoint,:);

                if any(isnan(startPoint))
                    warning('Intersect point missing')

                    %More start to top of aleurone, should be ok
                    startPoint = subscriptsToInterpolate(iPoint,:);

                    startAdjusted = 1;
                else
                    startAdjusted = 0; 
                end
            else
               % Endosperm, extend in from outer surface 
               startPoint = subscriptsToInterpolate(iPoint,:);
            end
        else
            % As above
            if pointIdentity
                %Aleurone, extend line in from interior surface
                startPoint = internalIntersectForSparse(iPoint-nPoints,:);

                if any(isnan(startPoint))
                    warning('Intersect point missing')

                    %Move start to top of aleurone, should be ok
                    startPoint = subscriptsForSparse(iPoint-nPoints,:);

                    startAdjusted = 1;
                else
                    startAdjusted = 0; 
                end
            else
               % Endosperm, extend in from outer surface 
               startPoint = subscriptsForSparse(iPoint-nPoints,:);
            end
        end

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

            if indexList(1) == tempIndexList(1) 

                if ~isempty(notEndospermID)
                   % Move start point back if possible

                    tempX = [fliplr(tempXNeg(2:(notEndospermID(1)))') tempX']';

                    tempY = [fliplr(tempYNeg(2:(notEndospermID(1)))') tempY']';

                    tempZ = [fliplr(tempZNeg(2:(notEndospermID(1)))') tempZ']';

                    indexList = [fliplr(tempIndexList(2:(notEndospermID(1)))') indexList']';

                    voxelDistances = [fliplr(tempDistancesNeg(2:(notEndospermID(1)))') voxelDistances']';
                else
                    warning('No precceding non-endosperm')
                end
            else
               error('Start indexes do not match') 
            end
        end

        lineIDs = grainVolumeAlignedCut(indexList);

        voxelSum = cumsum(voxelDistances);

        % Get endosperm voxels
        endospermPoints = find(lineIDs == ENDOSPERM_INDEX);

        if ~isempty(endospermPoints)

            % Endosperm should not be on first point
            if endospermPoints(1) == 1
               error('Endosperm point too early') 

               endospermPoints = [];
            end

            if pointIdentity

                % First point should always be al.
                if lineIDs(1) ~= ALEURONE_INDEX
                   error('First point is not al.') 

                   endospermPoints = [];
                end

                if ~startAdjusted
                    % All points up to first endosperm should be aleurone

                    if any(lineIDs(2:endospermPoints(1)-1) ~= ALEURONE_INDEX) 
                        warning('Not all al. before endosperm')

                        endospermPoints = [];
                    end
                else
                    % Will allow air at start in this case - sometimes breaks on germ

                    if any(lineIDs(2:endospermPoints(1)-1) ~= ALEURONE_INDEX & ...
                            lineIDs(2:endospermPoints(1)-1) ~= 0)
                        warning('Not all al. or air between surface and endosperm')

                        endospermPoints = [];
                    end
                end
            end

            % Remove endopserm after first 'real' break
            endospermDiff = diff(endospermPoints);

            endospermSteps = find(endospermDiff > 1);

            %[lineIDs(1:20) voxelSum(1:20)]

            % Fun note: This seems superficially simlar to switch debouncing
            % Note: By doing this correction for both aleurone thickness and
                % endosperm intensity, points on border could be counted in both if
                % 1 1 2 1 2 2 interlaced pattern. 
                % However, should be avoided after correction in aleurone tracer

            if ~isempty(endospermSteps)
                % Check if each step is less than half block size

                toInsert = [];
                for jStep = 1:length(endospermSteps)

                    testStart = endospermPoints(endospermSteps(jStep))+1;

                    testEnd = testStart - 1 + endospermDiff(endospermSteps(jStep))-1;

                    testInds = testStart:testEnd;

                    if sum(voxelDistances(testInds)) < blockThickness/2
                        %Remove if there is not a full block after
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

                         if all(lineIDs(testInds) == ENDOSPERM_INDEX)
                             % Insert as endosperm points
                             toInsert = [toInsert endospermPoints(endospermSteps(jStep))+1:...
                                 endospermPoints(endospermSteps(jStep)+1)-1];
                         else
                             % Remove following points
                             endospermPoints(endospermSteps(jStep)+1:end) = [];

                             warning('Part removed on following test')

                             break;
                         end
                    else
                        %Remove without test
                        endospermPoints(endospermSteps(jStep)+1:end) = [];

                        warning('Part removed on length test')

                        break
                    end
                end

                endospermPoints = sort([endospermPoints' toInsert]');
            end

            % Initalize lists with endosperm points
            indexList = indexList(endospermPoints);

            voxelDistances = (voxelDistances(endospermPoints));

            voxelSum = cumsum(voxelDistances);

            error('Indicate germ, al, air inclusions in depth profile')
            
            % Need at least one full block included
            if sum(voxelDistances) > blockThickness

                distanceIncluded = 0;

                tempProfile = zeros(1, numberOfBlocks)*NaN;

                % Step through blocks, averaging intensity
                for jBlock = 1:numberOfBlocks

                    % Take first set of blocks
                    blockSet = find(voxelSum < blockThickness & voxelDistances ~= 0);

                    if ~isempty(blockSet)
                        distanceSet = voxelDistances(blockSet);

                        %Check if final index in block is incomplete and there are following points
                        if voxelSum(blockSet(end)) < blockThickness & ...
                                blockSet(end) ~= length(indexList)
                            % If so, add next index after final to block
                            blockSet(end+1) = blockSet(end)+1;

                            %Add approrpaite distance, and calcualte amount remaining
                            distanceSet(end+1) = blockThickness - voxelSum(blockSet(end-1));

                            finalRemaining = voxelSum(blockSet(end)) - blockThickness;
                        else
                           %Block is perfect length
                           finalRemaining = 0; 
                        end

                        % Just include if thickness more than 3/4 of full blcok
                        if sum(distanceSet) > blockThickness*3/4;
                            % Take average
                            tempProfile(jBlock) = wmean( ...
                                double(greyVolumeAligned(indexList(blockSet))), distanceSet);

                            distanceIncluded = distanceIncluded + sum(distanceSet);
                        end

                        % Zero out remaining blocks and retake average
                        voxelDistances(blockSet) = 0;

                        if finalRemaining ~= 0
                            voxelDistances(blockSet(end)) = finalRemaining;
                        end

                        voxelSum = cumsum(voxelDistances);
                    end
                end

                % Save into correct array
                if iPoint <= nPoints
                    pointIntensityProfile(iPoint, :) = tempProfile;

                else
                    sparseIntensityProfile(iPoint - nPoints, :) = tempProfile;

                end


                if distanceIncluded < depthToCalculate
                   warning('Full depth not reached'); 
                end
            else
               warning('Less than one block') 
            end
        else
           warning('No endosperm, wrong normal?') 
        end
    else
       warning('No normal, skipped completely') 
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

plot(pointsUnwrapped((indexAbove),1), pointsUnwrapped((indexAbove),2), 'gd');

plot(pointsUnwrapped((indexBelow),1), pointsUnwrapped((indexBelow),2), 'cd');

%% Create nice image plotting

% Firstly create closed loop, similar problem as bee FOV
% Just based on angles to start with, can add TSP for more complicated shapes
sortedEdgeSubscripts = pointsUnwrapped(edgeIndexList, :); 
%createSortedLoopwithTSP(targetSubscripts(edgeIndexList, [1 3]));

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
% angles = atan2(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2));
        
% [~,  edgeIndexListSorted] = sort(angles);

% Shortest path solution, works really well!
[~, edgeIndexListSorted] = createSortedLoopwithTSP(double(sortedEdgeSubscripts));

% Add first point at end to close.
edgeIndexListSorted = [edgeIndexListSorted' edgeIndexListSorted(1)];

sortedEdgeSubscripts = sortedEdgeSubscripts(edgeIndexListSorted, :);

figure; hold on;

plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2))

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

plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2), 'b')

plot(offSetFullSubscripts(:,1), offSetFullSubscripts(:,2), 'rx')

plot(offSetSparseSubscripts(:,1), offSetSparseSubscripts(:,2), 'g.')

warning('Indicate crease edge')

%% Get voxels ID for aleurone and endosperm

IDInterpolant = scatteredInterpolant( double([offSetFullSubscripts(:,1)' offSetSparseSubscripts(:,1)']'), ...
    double([offSetFullSubscripts(:,2)' offSetSparseSubscripts(:,2)']'),...
    [interpolatedIdentity' sparseIdentity']',...
    'nearest','none');

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

%´Take borders

tempIm = IDImage;

% Just keep al to grow and get overlap
tempIm(IDImage == 1) = 0;

tempIm = imdilate(tempIm, strel('disk',1));

borderInds = find(IDImage == 1 & tempIm);

[borderX, borderY] = ind2sub([xRange , yRange], borderInds);

hold on

plot(borderX, borderY, 'g.')

warning('Get full border from image as for endopserm intersect')

%% Calculate thickness and intesntiy image
%figure; subplot(1,2,1); hist(thicknessByPoint,100)
%subplot(1,2,2); hist(averageIntensityByPoint,100)

warning('add voxel size to thickness scaling')
valuesToUse = find(~isnan(thicknessByPoint));

valuesToSparse = find(~isnan(thicknessForSparse));

thicknessInterpolant = scatteredInterpolant( double([offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']'), ...
    double([offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']'),...
    [thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']',...
    'linear','none');

tempImage = zeros(xRange , yRange);

tempImage(inAleurone) = thicknessInterpolant(XPointsInAleurone, YPointsInAleurone);

thicknessImage = zeros(xRange , yRange, 3);

thicknessImage(:,:,1) = tempImage;

warning('Colour range setting is not automated')

[max(thicknessByPoint) max(thicknessForSparse)]

%thicknessImage = (thicknessImage-2*VOXEL_SIZE)/(18*VOXEL_SIZE-2*VOXEL_SIZE);
warning('Change thicknes scaling')
thicknessImage = (thicknessImage-2)/(18-2);

endoImage = zeros(xRange , yRange);

endoImage(inEndosperm) = 0.4;

%thicknessImage(:,:,2) = endoImage;

cols = zeros(100,3);
cols(1:100,1) = (1:100)/100;

figure; 
imshow(permute(thicknessImage, [2 1 3]))
colormap(cols)
title('Thickness'); hcb = colorbar; set(hcb,'Ticks', [0 0.5 1], 'TickLabels', {'8','40','72'})

hold on
plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2), 'w.')
plot(borderX, borderY, 'g.')

% Do for intensity.
intensityInterpolant = scatteredInterpolant( double([offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']'), ...
    double([offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']'),...
    [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']',...
    'linear','none');

tempImage = zeros(xRange , yRange);

tempImage(inAleurone) = intensityInterpolant(XPointsInAleurone, YPointsInAleurone);

intensityImage = zeros(xRange , yRange, 3);

intensityImage(:,:,3) = tempImage;

intensityImage = tempImage;

[max(averageIntensityByPoint) max(averageIntensityForSparse)]
[min(averageIntensityByPoint) min(averageIntensityForSparse)]

intensityImage = (intensityImage-50)/(250-50);

%intensityImage(:,:,2) = endoImage;

cols = zeros(100,3);
cols(1:100,3) = (1:100)/100;

figure;
imshow(permute(intensityImage, [2 1 3]))
colormap(cols)
title('Intensity'); hcb = colorbar; set(hcb,'Ticks', [0 0.5 1], 'TickLabels', {'50','150','250'})

hold on
plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2), 'w.')
plot(borderX, borderY, 'g.')

% Make overlay.
overlayImage = thicknessImage + intensityImage;

overlayImage(:,:,2) = overlayImage(:,:,2)/2;

figure;
imshow(permute(overlayImage, [2 1 3]));
title('Overlay');

hold on
plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2), 'w.')
plot(borderX, borderY, 'g.')


figure; 
plot([thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']', ...
    [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']', '.');
rValue = corr([thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)']', ...
    [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)']');
title(sprintf('%.2f', rValue))
%% Calculate and display intensity depth maps
figure; imshow(intensityImage')

error('Create single 3D interpolant and block on inclusions')

figure;

[max(pointIntensityProfile(:)) max(sparseIntensityProfile(:))]
[min(pointIntensityProfile(:)) min(sparseIntensityProfile(:))]

for iBlock = 1:numberOfBlocks
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock)));
    
    intensityProfileInterpolant = scatteredInterpolant( double([offSetFullSubscripts(valuesToUse,1)' offSetSparseSubscripts(valuesToSparse,1)']'), ...
        double([offSetFullSubscripts(valuesToUse,2)' offSetSparseSubscripts(valuesToSparse,2)']'),...
        [pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)']',...
        'linear','none');
    
    tempImage = zeros(xRange , yRange);
    
    tempImage(inMap) = intensityProfileInterpolant(XPointsIn, YPointsIn);
    
    tempImage = (tempImage-0)/(100-0);
    
    subplot(2,5,iBlock)
    imshow(tempImage'); hold on;
    
    plot(sortedEdgeSubscripts(:,1), sortedEdgeSubscripts(:,2), 'w.')
    plot(borderX, borderY, 'g.')
end
    
%% Look at intensity correlation between aleurone and depth layers
figure;

for iBlock = 1:numberOfBlocks
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)) & ~isnan(thicknessByPoint));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock)) & ~isnan(thicknessForSparse));
    
    alIntensity = [averageIntensityByPoint(valuesToUse)' averageIntensityForSparse(valuesToSparse)'];
    
    blockIntensity = [pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)'];
    
    subplot(2,5,iBlock)
    plot(alIntensity, blockIntensity, '.');
    
    rValue = corr(alIntensity, blockIntensity);
    title(sprintf('%.2f', rValue))
    
end

figure;

for iBlock = 1:numberOfBlocks
    valuesToUse = find(~isnan(pointIntensityProfile(:,iBlock)) & ~isnan(thicknessByPoint));

    valuesToSparse = find(~isnan(sparseIntensityProfile(:,iBlock)) & ~isnan(thicknessForSparse));
    
    alThickness = [thicknessByPoint(valuesToUse)' thicknessForSparse(valuesToSparse)'];
    
    blockIntensity = [pointIntensityProfile(valuesToUse,iBlock)' sparseIntensityProfile(valuesToSparse,iBlock)'];
    
    subplot(2,5,iBlock)
    plot(alIntensity, blockIntensity, '.');
    
    rValue = corr(alIntensity, blockIntensity);
    title(sprintf('%.2f', rValue))
    
end

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
