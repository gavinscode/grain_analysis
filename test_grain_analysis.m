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

tempIndexList = sub2ind(volumeSize, curveSubscriptArray(:,1), curveSubscriptArray(:,2), curveSubscriptArray(:,3));

% Code curve into grain exterior as 2 to distinguish it later.
grainExterior = uint8(grainExterior);
grainExterior(tempIndexList(grainExterior(tempIndexList) == 1)) = 2;

figure; imshow(sum(grainExterior(:,:,1623),3)/2)

clear centreCurveVolume, 

%% Get exterior surface of aleurone.
aleuroneExterior = (grainExterior == 1) & (grainVolumeAligned == ALEURONE_INDEX);

% Take largest connected region of surface.
tempCC = bwconncomp(aleuroneExterior, 26);

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
    aleuroneExterior(tempStats(iRegion).PixelIdxList) = 0;
end



% Get edge of aleurone surface.
aleuroneEdge = imdilate(grainExterior & (grainVolumeAligned == ENDOSPERM_INDEX | ...
    grainVolumeAligned == GERM_INDEX | grainExterior == 2), STREL_18_CONNECTED);

aleuroneEdge = aleuroneEdge & aleuroneExterior;

% Take largest connected region of edge.
tempCC = bwconncomp(aleuroneEdge, 26);

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
    aleuroneEdge(tempStats(iRegion).PixelIdxList) = 0;
end


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
% plot3(aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), ...
%     aleuroneSurfaceSubscriptArray(:,3), 'b.')
plot3(aleuroneEdgeSubscriptArray(:,1), aleuroneEdgeSubscriptArray(:,2), ...
    aleuroneEdgeSubscriptArray(:,3), 'r.')

[sum(aleuroneExterior(aleuroneSurfaceIndexList)) sum(aleuroneExterior(aleuroneEdgeIndexList))]



% Get endosperm and germ surfaces for plotting.
endospermSurfaceIndexList = find(grainExterior & (grainVolumeAligned == ENDOSPERM_INDEX));

nIndex = length(endospermSurfaceIndexList); endospermSurfaceSubscriptArray = zeros(nIndex, 3);

[endospermSurfaceSubscriptArray(:,1), endospermSurfaceSubscriptArray(:,2), endospermSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, endospermSurfaceIndexList);

germSurfaceIndexList = find(grainExterior & (grainVolumeAligned == GERM_INDEX));

nIndex = length(germSurfaceIndexList); germSurfaceSubscriptArray = zeros(nIndex, 3);

[germSurfaceSubscriptArray(:,1), germSurfaceSubscriptArray(:,2), germSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, germSurfaceIndexList);

clear aleuroneEdge, clear grainExterior
%% Allocate points equally across surface and border, as in bee
edgePointsToChoose = 1:length(aleuroneEdgeIndexList);

surfacePointsToChoose = 1:length(aleuroneSurfaceIndexList);

edgePointsChoosen = zeros(length(aleuroneEdgeIndexList),1);

surfacePointsChoosen = zeros(length(aleuroneSurfaceIndexList),1);

% for 300, 500 gives 12 9 points - ~20 minutes, but this does not unwrap nicely
% for 20, 50, gives 343 1314 points - ~27 hours
% for 50, 100, gives 121 309 points - ~7 hours
% for 50, 75, gives 121 565 points - ~11 hours - go with this!
% for 30, 75, gives 216 562 points - ~13 hours 

% For test, seems good to get dense on border and moderate on surface

edgeDistance = 25; 

surfaceDistance = 500;

% Test edge points first. Select from top of germ (max Y)
while ~isempty(edgePointsToChoose)
    [~, ind] = max(aleuroneEdgeSubscriptArray(edgePointsToChoose,3));
    
    edgePointsChoosen(edgePointsToChoose(ind)) = 1;
    
    pointChoosen = aleuroneEdgeSubscriptArray(edgePointsToChoose(ind),:);
    
    % Remove edge and surface points nearby.
    indsToRemove = find( sqrt( (aleuroneEdgeSubscriptArray(edgePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (aleuroneEdgeSubscriptArray(edgePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (aleuroneEdgeSubscriptArray(edgePointsToChoose,3) - pointChoosen(3)).^2) < edgeDistance);
    
    edgePointsToChoose(indsToRemove) = [];
    
    indsToRemove = find( sqrt( (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,3) - pointChoosen(3)).^2) < surfaceDistance);
    
    surfacePointsToChoose(indsToRemove) = [];
end

% Select surface points from remaining, again going down Z
while ~isempty(surfacePointsToChoose)
    [~, ind] = max(aleuroneSurfaceSubscriptArray(surfacePointsToChoose,3));
    
    surfacePointsChoosen(surfacePointsToChoose(ind)) = 1;
    
    pointChoosen = aleuroneSurfaceSubscriptArray(surfacePointsToChoose(ind),:);
    
    % Remove surface points nearby.  
    indsToRemove = find( sqrt( (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,1) - pointChoosen(1)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,2) - pointChoosen(2)).^2 + ...
        (aleuroneSurfaceSubscriptArray(surfacePointsToChoose,3) - pointChoosen(3)).^2) < surfaceDistance);
    
    surfacePointsToChoose(indsToRemove) = [];
end

edgePointsChoosen = find(edgePointsChoosen); 

surfacePointsChoosen = find(surfacePointsChoosen); 

[length(edgePointsChoosen) length(surfacePointsChoosen)]

figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
plot3(aleuroneSurfaceSubscriptArray(surfacePointsChoosen,1), aleuroneSurfaceSubscriptArray(surfacePointsChoosen,2), ...
    aleuroneSurfaceSubscriptArray(surfacePointsChoosen,3), 'b.')
plot3(aleuroneEdgeSubscriptArray(edgePointsChoosen,1), aleuroneEdgeSubscriptArray(edgePointsChoosen,2), ...
    aleuroneEdgeSubscriptArray(edgePointsChoosen,3), 'r.')
%% Calulate distance between points

%%% Test how is geodesic connectivity calculated - its 26
%temp = zeros(3,3,3,'logical');
%temp(1, :, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 1) = 1; temp(3, 3, 1) = 1;
%temp(1, 1, 1) = 1; temp(2, 2, 2) = 1; temp(3, 3, 3) = 1;
%bwdistgeodesic(temp, 1,'quasi-euclidean')

subscriptsToInterpolate = [aleuroneSurfaceSubscriptArray(surfacePointsChoosen,:)'...
    aleuroneEdgeSubscriptArray(edgePointsChoosen,:)']';

indsToInterpolate = sub2ind(volumeSize, subscriptsToInterpolate(:,1), subscriptsToInterpolate(:,2),...
    subscriptsToInterpolate(:,3));

nPoints = length(indsToInterpolate);

distanceMatrix = zeros(nPoints, nPoints);

%%% Should be able to parallelize this.

%%% Load in grey scale volume here for average calculations. 

% Loop through geodesic distance calculations for each point than put into matrix.
for iPoint = 1; %:nPoints %
    tic
    dMap = bwdistgeodesic(aleuroneExterior, indsToInterpolate(iPoint),'quasi-euclidean');
    toc

    %tempInd = sub2ind(volumeSize, 378, 175, 1651);
    %dMap = bwdistgeodesic(aleuroneExterior, tempInd,'quasi-euclidean');
    
    warning('Get normal, average intensity, thickness')
    %%% Check normal code carefully...
    
    % Pause to let matlab free memory (?)
    pause(0.1)
    
    distanceMatrix(iPoint, :) = dMap(indsToInterpolate);
end

% save(sprintf('/Users/gavintaylor/Documents/Matlab/Temp_data/%s_%i_%i', 'distanceMatrix', edgeDistance, surfaceDistance), 'distanceMatrix');
load('/Users/gavintaylor/Documents/Matlab/Temp_data/distanceMatrix_50_75.mat')

% Save to prevent overwrite.
% distanceMatrixSaver = distanceMatrix;
% pPoint = iPoint;

% Debug plot - show colour maps. Tests well on 7
figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
cols = round(dMap(aleuroneSurfaceIndexList)/max(distanceMatrix(:))*99) + 1;
cols(isinf(cols)) = 110;
fscatter3(aleuroneSurfaceSubscriptArray(:,1),aleuroneSurfaceSubscriptArray(:,2),...
    aleuroneSurfaceSubscriptArray(:,3),cols,jet(100));
% plot3(subscriptsToInterpolate(pPoint,1), subscriptsToInterpolate(pPoint,2), ...
%     subscriptsToInterpolate(pPoint,3),'mo', 'markersize',20);
% plot3(378,175,1651,'mo', 'markersize',20);
plot3(subscriptsToInterpolate(iPoint,1), subscriptsToInterpolate(iPoint,2), ...
    subscriptsToInterpolate(iPoint,3),'mo', 'markersize',20);


% Debug plot - show series of slices 
% dMap(isnan(dMap)) = 0;
% for iTemp = 1615:1625
%     figure; imshow(sum(dMap(:,:,iTemp:iTemp+1),3))
% end
% 
% %Debug plot - show thick overlay figure
% figure; imshow(sum(dMap(:,:,1623:1624),3))
% 
% % Debug plot - show local series of points
% inds = find(aleuroneSurfaceSubscriptArray(:,3) > 1615 & aleuroneSurfaceSubscriptArray(:,3) < 1625);
% figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
% plot3(aleuroneSurfaceSubscriptArray(inds,1),aleuroneSurfaceSubscriptArray(inds,2),...
%      aleuroneSurfaceSubscriptArray(inds,3),'b.');
%  
% % %Debug plot - show series of coloured points
% inds = find(aleuroneSurfaceSubscriptArray(:,3) > 1600 & aleuroneSurfaceSubscriptArray(:,3) < 1700 & ...
%     aleuroneSurfaceSubscriptArray(:,1) < 380);
% inds = find(aleuroneSurfaceSubscriptArray(:,3) < 200);
% cols = round(dMap(aleuroneSurfaceIndexList(inds))/max(dMap(aleuroneSurfaceIndexList(inds)))*99) + 1;
% cols(isinf(cols)) = 110; cols(cols > 100) = 100;
% figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
% fscatter3(aleuroneSurfaceSubscriptArray(inds,1),aleuroneSurfaceSubscriptArray(inds,2),...
%      aleuroneSurfaceSubscriptArray(inds,3),cols,jet(100));

%Just test lines between points, can do on any afterwards
% figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
% plot3(aleuroneSurfaceSubscriptArray(surfacePointsChoosen,1), aleuroneSurfaceSubscriptArray(surfacePointsChoosen,2), ...
%     aleuroneSurfaceSubscriptArray(surfacePointsChoosen,3), 'b.')
% plot3(aleuroneEdgeSubscriptArray(edgePointsChoosen,1), aleuroneEdgeSubscriptArray(edgePointsChoosen,2), ...
%     aleuroneEdgeSubscriptArray(edgePointsChoosen,3), 'r.')
% plot3(aleuroneEdgeSubscriptArray(:,1), aleuroneEdgeSubscriptArray(:,2), ...
%     aleuroneEdgeSubscriptArray(:,3), 'g.')
% 
% colMap = jet(100);
% maxDist = max(distanceMatrix(:));
% 
% pointToTest = 23; % prob on 20 to 7
% 
% for iPoint = 1:nPoints
%     col = colMap(round(distanceMatrix(pointToTest,iPoint)/maxDist*99)+1,:);
%     line([subscriptsToInterpolate(pointToTest,1) subscriptsToInterpolate(iPoint,1)], ...
%         [subscriptsToInterpolate(pointToTest,2) subscriptsToInterpolate(iPoint,2)], ...
%         [subscriptsToInterpolate(pointToTest,3) subscriptsToInterpolate(iPoint,3)], 'color', col);
% end

% plot3(subscriptsToInterpolate(pPoint,1), subscriptsToInterpolate(pPoint,2), ...
%     subscriptsToInterpolate(pPoint,3),'mo', 'markersize',20);

clear aleuroneExterior, clear dMap

%% Calculate unwrapping - based on tutorial on Numerical tours
% https://nbviewer.jupyter.org/github/gpeyre/numerical-tours/blob/master/matlab/meshdeform_3_flattening.ipynb

% Check for points with inf distance.
if any(isinf(distanceMatrix(:))) || any(isnan(distanceMatrix(:)))
   error('Can not have NaN or Inf in distance matrix')
end

% Enforce symmetry (should be ok...)
distanceMatrixTemp = (distanceMatrix + distanceMatrix')/2;

% Compute centered matrix.
J = eye(nPoints) - ones(nPoints)/nPoints;
W = -J*(distanceMatrixTemp.^2)*J;

% Diagonalize centred matrix.
[U,S] = eig(W);
S = diag(S);
[S,I] = sort(S,'descend'); 
U = U(:,I);

figure; plot(S,'x-')

% Map is determined by two largest values
%%% Take three largest if considering manifold
pointsUnwrapped = (U(:,1:2)' .* repmat(sqrt(S(1:2)), [1 nPoints]))';

pointsUnwrapped = pointsUnwrapped/2;

%%% Should check how well diemsions are preserved during unwrapping.
% D2 will go to Z, D1 will go to X, Y is set in middle

%%% Need to majorly reshape new points to preserve dimensions. Done manually for now
% Flip Y and add minimum
pointsUnwrapped(:,2) = -pointsUnwrapped(:,2); 
pointsUnwrapped(:,2) = pointsUnwrapped(:,2) - min(pointsUnwrapped(:,2));

%Flip X, center
pointsUnwrapped(:,1) = -pointsUnwrapped(:,1);
pointsUnwrapped(:,1) = pointsUnwrapped(:,1)-mean(pointsUnwrapped(:,1))+volumeSize(1)/2;

targetSubscripts = [pointsUnwrapped(:,1) ones(size(pointsUnwrapped,1),1)*volumeSize(2)/2, pointsUnwrapped(:,2)];

figure; 
subplot(1, 2, 1); hold on; axis equal; set(gca, 'Clipping', 'off'); axis off

plot3(endospermSurfaceSubscriptArray(:,1), endospermSurfaceSubscriptArray(:,2), ...
    endospermSurfaceSubscriptArray(:,3), 'g.')

plot3(germSurfaceSubscriptArray(:,1), germSurfaceSubscriptArray(:,2), ...
    germSurfaceSubscriptArray(:,3), 'b.')

plot3(aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), ...
    aleuroneSurfaceSubscriptArray(:,3), 'y')

plot3(aleuroneSurfaceSubscriptArray(surfacePointsChoosen,1), aleuroneSurfaceSubscriptArray(surfacePointsChoosen,2), ...
    aleuroneSurfaceSubscriptArray(surfacePointsChoosen,3), 'bo')

plot3(aleuroneEdgeSubscriptArray(edgePointsChoosen,1), aleuroneEdgeSubscriptArray(edgePointsChoosen,2), ...
    aleuroneEdgeSubscriptArray(edgePointsChoosen,3), 'ro')

% for i = 1:length(edgePointsChoosen)
%     text(aleuroneEdgeSubscriptArray(edgePointsChoosen(i),1), aleuroneEdgeSubscriptArray(edgePointsChoosen(i),2), ...
%     aleuroneEdgeSubscriptArray(edgePointsChoosen(i),3), sprintf('%i', i ));
% end

subplot(1, 2, 2); hold on; axis equal; set(gca, 'Clipping', 'off'); axis off
plot3(targetSubscripts(1:length(surfacePointsChoosen),1), targetSubscripts(1:length(surfacePointsChoosen),2), ...
    targetSubscripts(1:length(surfacePointsChoosen),3), 'b.')

plot3(targetSubscripts(length(surfacePointsChoosen)+1:end,1), targetSubscripts(length(surfacePointsChoosen)+1:end,2), ...
    targetSubscripts(length(surfacePointsChoosen)+1:end,3), 'r.')

warning('Add interpolation between points, check gnat')

%%% Get distance bar by calculating distance between pairs of points in 2D
%%% and then comparing to geodesic, take average as scale (variance?)
warning('Add scale bar');

% for i = length(surfacePointsChoosen)+1:size(targetSubscripts,1)
%     text(targetSubscripts(i,1), targetSubscripts(i,2), targetSubscripts(i,3), sprintf('%i', i - length(surfacePointsChoosen)));
% end


%% Now use point_registration and bspline_transform to shift volumes.
%%%% May wish to draw lines in from points (based on surface normal) to get
%%% interior manifold. These points are then just measured distance beneath
%%% others.

%%% Other trick, try to do with interpolation. 
% Get dense points on outer and inner surface. 
% Have coordinates along linkage path.

options.MaxRef=5;
options.Verbose=true;

% Simply scaling down points and volume size speeds regularization up a lot...
%%% Creates same results on reduced size image. May be better way to scale up?
sF = 4;
[O_trans, Spacing] = point_registration(round(volumeSize([1 2 3])/sF), subscriptsToInterpolate(:, [1 2 3])/sF, ...
    targetSubscripts(:, [1 2 3])/sF, options);

%%% Diffomorphic is require to prevent problems, but takes forever on full volume.
    % Try with one it of point placing (coase spacing). A bit faster...
%%% Ok Regularize does not work, probably means b-spline is not goind to work.
   % Doing manifold (with 3 dimensions from parameterization) may help
   % enforce order on spline.
 
% tic
% [O_trans,Spacing]=MakeDiffeomorphic(O_trans,Spacing,round(volumeSize([1 2 3])/sF));
% toc

%%% Target points match exactly, other points are generally ok and lie on plane
    %%% Although the border is quite wavy 
    %%% and some points cross between flanges near germ.
       %%% sorting borders points so their order doesn't cross back and forth could help.
       %%% Cross over occurs where misplaced points are...
%%% Dense point sets helps a lot, but internal manifold will probably not assist with this.
tic
flatEdgeSubscriptArray = bspline_trans_points_double(O_trans, Spacing, ...
    aleuroneEdgeSubscriptArray/sF)*sF;

flatSurfaceSubscriptArray = bspline_trans_points_double(O_trans, Spacing, ...
    aleuroneSurfaceSubscriptArray(1:100:end,:)/sF)*sF;
    %
toc

%For Z 1000....
% Use directly 
    % Stretches along crease
% Flip X and Y on O trans (and volume size) 
    % Opens on one side but rotates 90 (has dent back opposite crease)
% Flip X and Y on subscripts input (and volume size), but not on O trans.
    % Opens on one side (possibly flipped, and dent back opposite crease)
% Flip X and Y on subscripts input (and volume size), and on O trans.
    % Rotated 90 and stretech along grain
% Flip X and Y on subscripts input (and volume size), but not on O trans, and volume permuted.
    % As direct...

%%% Many strange artifacts using spline on image, diffeomrophic regularization may help
%%% Placing internal manifold may also help.
% tic
% flatGrain = uint8( bspline_transform(O_trans(:, :, :, [1 2]), ... 
%    single( grainVolumeAligned(:,:,1400)), Spacing([1 2])));
% % flatGrain = permute(uint8( bspline_transform(O_trans(:, :, :, :), ... 
% %    single( permute(grainVolumeAligned, [2 1 3])), Spacing)), [2 1 3]);
% toc

% slice = 1400;
% figure; subplot(1,2,1); imshow(grainVolumeAligned(:,:,slice)*50); hold on;
% indsToPlot = find(aleuroneSurfaceSubscriptArray(:,3) == slice);
% plot(aleuroneSurfaceSubscriptArray(indsToPlot,2), aleuroneSurfaceSubscriptArray(indsToPlot,1),'r.')
% subplot(1,2,2); imshow(flatGrain*50)
%imshow(flatGrain(:,:,slice)*50)

% Test unwrapping of point cloud
figure; hold on ; axis equal; set(gca, 'Clipping', 'off')
plot3(targetSubscripts(:,1), targetSubscripts(:,2), ...
    targetSubscripts(:,3), 'rx')

plot3(aleuroneEdgeSubscriptArray(:,1), aleuroneEdgeSubscriptArray(:,2), ...
    aleuroneEdgeSubscriptArray(:,3), 'r.')

%plot3(flatSubscriptArray(1:100:end,1), flatSubscriptArray(1:100:end,2), ...
%    flatSubscriptArray(1:100:end,3), 'g.')

plot3(flatEdgeSubscriptArray(:,1), flatEdgeSubscriptArray(:,2), ...
    flatEdgeSubscriptArray(:,3), 'mo')

plot3(flatSurfaceSubscriptArray(:,1), flatSurfaceSubscriptArray(:,2), ...
    flatSurfaceSubscriptArray(:,3), 'bo')