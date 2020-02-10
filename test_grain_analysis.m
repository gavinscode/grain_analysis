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
%%% This can be solved by doing sequentially for small sequences and then recombining.
%%% Do say, for 1-50 and 50-100 then combine. Should test 25-75 produces same result in overlap.
warning('Not using all points in nrbloft')
tic
curveNrb = nrbloft_nan(xArray(:,1:10:end), yArray(:,1:10:end), zArray(:,1:10:end), 2);
toc

%% Put interpolated values into curve volume.
%zToInterp = zBottomOfLoop-zTopOfLoop+1;
zToInterp = max(max(zArray(:,1:10:end)))-zTopOfLoop+1;

%yToInterp = max(loopSubscriptArray(:,2));
yToInterp = max(max(yArray(:,1:10:end)));

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
centreCurveVolume = centreCurveVolume - imerode(centreCurveVolume, STREL_18_CONNECTED);

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

curveIndexList = find(centreCurveVolume);

nIndex = length(curveIndexList); curveSubscriptArray = zeros(nIndex, 3);

[curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3)] = ...
    ind2sub(smallVolumeSize, curveIndexList);

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
% Correct subscripts back to original coordinates
curveSubscriptArray(:,1) = curveSubscriptArray(:,1) + xBoundsNew(1) - 1;

curveSubscriptArray(:,2) = curveSubscriptArray(:,2) + yBoundsNew(1) - 1;

curveSubscriptArray(:,3) = curveSubscriptArray(:,3) + zBoundsNew(1) - 1;

tempIndexList = sub2ind(volumeSize, curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3));

grainExterior(tempIndexList) = 0;

%clear centreCurveVolume

%% Get exterior of aleurone.
aleuroneExterior = grainExterior & (grainVolumeAligned == ALEURONE_INDEX);

% Get edge of aleurone exterior.
aleuroneEdge = imdilate(grainExterior & (grainVolumeAligned == ENDOSPERM_INDEX | ...
    grainVolumeAligned == GERM_INDEX), STREL_18_CONNECTED);

aleuroneEdge = aleuroneEdge & aleuroneExterior;

aleuroneExterior = aleuroneExterior - aleuroneEdge;

% Get index list and test plot both.
aleuroneSurfaceIndexList = find(aleuroneExterior);

nIndex = length(aleuroneSurfaceIndexList); aleuroneSurfaceSubscriptArray = zeros(nIndex, 3);

[aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), aleuroneSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, aleuroneSurfaceIndexList);

aleuroneEdgeIndexList = find(aleuroneEdge);

nIndex = length(aleuroneEdgeIndexList); aleuroneEdgeSubscriptArray = zeros(nIndex, 3);

[aleuroneEdgeSubscriptArray(:,1), aleuroneEdgeSubscriptArray(:,2), aleuroneEdgeSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, aleuroneEdgeIndexList);

figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
plot3(aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), ...
    aleuroneSurfaceSubscriptArray(:,3), 'b.')
plot3(aleuroneEdgeSubscriptArray(:,1), aleuroneEdgeSubscriptArray(:,2), ...
    aleuroneEdgeSubscriptArray(:,3), 'r.')

clear aleuroneExterior, clear aleuroneEdge,  

%% Allocate points equally across surface and border, as in bee
edgePointsToChoose = 1:length(aleuroneEdgeIndexList);

surfacePointsToChoose = 1:length(aleuroneSurfaceIndexList);

edgePointsChoosen = zeros(length(aleuroneEdgeIndexList),1);

surfacePointsChoosen = zeros(length(aleuroneSurfaceIndexList),1);

edgeDistance = 50;

surfaceDistance = 100;

% Test edge points first. Select from top of germ (max Y)
while ~isempty(edgePointsToChoose)
    [~, ind] = max(aleuroneEdgeSubscriptArray(edgePointsToChoose,3));
    
    edgePointsChoosen(ind) = 1;
    
    pointChoosen = aleuroneEdgeSubscriptArray(edgePointsToChoose(ind),:);
    
    % Remove edge and surface points nearby.
    indsToRemove = find( sqrt( (aleuroneEdgeSubscriptArray(edgePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (aleuroneEdgeSubscriptArray(edgePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (aleuroneEdgeSubscriptArray(edgePointsToChoose,3) - pointChoosen(3)).^2) < edgeDistance);
    
    edgePointsToChoose(indsToRemove) = [];
    
    indsToRemove = find( sqrt( (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,3) - pointChoosen(3)).^2) < edgeDistance);
    
    surfacePointsToChoose(indsToRemove) = [];
end

error('Find out why surface points not allocated on bottom of grain')

% Select surface points from remaining, now go across X
while ~isempty(surfacePointsToChoose)
    [~, ind] = max(aleuroneSurfaceSubscriptArray(surfacePointsToChoose,1));
    
    surfacePointsChoosen(ind) = 1;
    
    pointChoosen = aleuroneSurfaceSubscriptArray(surfacePointsToChoose(ind),:);
    
    % Remove surface points nearby.  
    indsToRemove = find( sqrt( (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,3) - pointChoosen(3)).^2) < surfaceDistance);
    
    surfacePointsToChoose(indsToRemove) = [];
end

edgePointsChoosen = find(edgePointsChoosen);

surfacePointsChoosen = find(surfacePointsChoosen);

figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
plot3(aleuroneSurfaceSubscriptArray(surfacePointsChoosen,1), aleuroneSurfaceSubscriptArray(surfacePointsChoosen,2), ...
    aleuroneSurfaceSubscriptArray(surfacePointsChoosen,3), 'b.')
plot3(aleuroneEdgeSubscriptArray(edgePointsChoosen,1), aleuroneEdgeSubscriptArray(edgePointsChoosen,2), ...
    aleuroneEdgeSubscriptArray(edgePointsChoosen,3), 'r.')

%% Calulate distance points

%%% Test how is geodesic connectivity calculated - its 26
%temp = zeros(3,3,3,'logical');
%temp(1, :, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 1) = 1; temp(3, 3, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 2) = 1; temp(3, 3, 3) = 1;
%bwdistgeodesic(temp, 1,'quasi-euclidean')

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