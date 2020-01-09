% Test load and plot of grain stack.
clc; clear; close all

% Other: Om_1_7_test
targetDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Grain LU/Data/OB6_test/Labels';

% Loading tiff stack takes a while.
grainVolume = loadtiffstack(targetDirectory, 1);

grainVolume = grainVolume(:,:,1:1500);

volumeSize = size(grainVolume);

aleuroneIndex = 1; % Material index.



%% Take mesh around entire grain then limit to aleurone exterior.

% Reduce patch size to save mesh construction time - keep as large as possible.
% If this value is below ~0.05, holes may appear in mesh if grain touches volume border.
resizeRatio = 0.2;

% Do a strong morphological open to seperate along crease and remove floaters.
%smallGrainVolume = imopen(grainVolume, strel('disk', 5));
smallGrainVolume = grainVolume;

% Volume should not be binarized before this point.
smallGrainVolume = imresize3(uint8(smallGrainVolume), resizeRatio);

%%% Should test for presence of multiple objects/floaters here.

% Create isosurface on downsampled volume.
fullGrainSurface = isosurface(smallGrainVolume, 0.5);

%voxelizedFullGrainSurface = inpolyhedron(fullGrainSurface, 1:voxelSize(1), ... 
%   1:voxelSize(2), 1:voxelSize(3));

% Simplify patch to speed processing time later. 
% This and resize ratio both need tuning to result in < 5000 vertices.
patchReductionRatio = 0.02;

fullGrainSurface = reducepatch(fullGrainSurface,patchReductionRatio);

length(fullGrainSurface.vertices)

% Englarge surface back to original size.
fullGrainSurface.vertices = round(fullGrainSurface.vertices/resizeRatio);

% Note that surface vertices are flipped relative to image subscripts.
fullGrainSurface.vertices(:,[1 2]) = fullGrainSurface.vertices(:,[2 1]);

% Was testing mesh voxilization for debugging
%voxelizedSmallGrainSurface = inpolyhedron(fullGrainSurface, 1:volumeSize(1), ... 
%   1:volumeSize(2), 1:volumeSize(3));

% Truncate to lower bound of volume.
fullGrainSurface.vertices(fullGrainSurface.vertices < 1) = 1;
        
% Truncate to upper bound.
fullGrainSurface.vertices(fullGrainSurface.vertices(:,1) > volumeSize(1),1) = volumeSize(1);

fullGrainSurface.vertices(fullGrainSurface.vertices(:,2) > volumeSize(2),2) = volumeSize(2);

fullGrainSurface.vertices(fullGrainSurface.vertices(:,3) > volumeSize(3),3) = volumeSize(3);

nVertex = size(fullGrainSurface.vertices,1);

nFace = size(fullGrainSurface.faces,1);

% Get long axis of grain with PCA.
grainAxisArray = pca(fullGrainSurface.vertices);

grainLongAxis = grainAxisArray(:,1)';

grainCreaseAxis = grainAxisArray(:,3)';

grainCenter = mean(fullGrainSurface.vertices);

% Plost test figures;
figure; axis equal; hold on; axis off; set(gca, 'Clipping', 'off')

patch(fullGrainSurface);

plot3(fullGrainSurface.vertices(:,1), fullGrainSurface.vertices(:,2), fullGrainSurface.vertices(:,3), 'ro')

line(grainCenter(1)+[0 2000]*grainLongAxis(1), grainCenter(2)+[0 2000]*grainLongAxis(2), ...
    grainCenter(3)+[0 2000]*grainLongAxis(3), 'color' ,'b')

line(grainCenter(1)+[0 1000]*grainAxisArray(1,2), grainCenter(2)+[0 1000]*grainAxisArray(2,2), ...
    grainCenter(3)+[0 1000]*grainAxisArray(3,2), 'color' ,'g')

line(grainCenter(1)+[0 1000]*grainCreaseAxis(1), grainCenter(2)+[0 1000]*grainCreaseAxis(2), ...
    grainCenter(3)+[0 1000]*grainCreaseAxis(3), 'color' ,'r')

title('Check there are no holes in the surface')



%% Get voxels on exterior of orignal grain by overlap to exterior.
%%% ADD FUNCTION for this and replace in loop.
grainExterior = ~grainVolume;

grainExterior = imdilate(grainExterior, strel('disk', 1));

grainExterior = grainExterior & grainVolume;

% Get surface indices, then convert to subscripts.
voxelIndexList = find(grainExterior);

nIndex = length(voxelIndexList);

grainSurfaceSubscriptArray = zeros(nIndex, 3);

[grainSurfaceSubscriptArray(:,1), grainSurfaceSubscriptArray(:,2), grainSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, voxelIndexList);



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
    
    if aleuroneIndex ~= voxelMaterialIndex
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
transform2Vertical = matrix2rotatevectors([0, 0, 1], grainLongAxis);

aleuroneSurfaceRotated = aleuroneSurface;

aleuroneSurfaceRotated.vertices = aleuroneSurfaceRotated.vertices - grainCenter;

aleuroneSurfaceRotated.vertices = aleuroneSurfaceRotated.vertices*transform2Vertical;

transform2Up = matrix2rotatevectors([0, 1, 0], aleuroneCreaseAxis*transform2Vertical);

aleuroneSurfaceRotated.vertices = aleuroneSurfaceRotated.vertices*transform2Up;

% Do for Aleurone surface.
% voxelIndexList = find(grainExterior & grainVolume == 1);
% 
% nIndex = length(voxelIndexList);
% 
% aleuroneSurfaceSubscriptArray = zeros(nIndex, 3);
% 
% [aleuroneSurfaceSubscriptArray(:,1), aleuroneSurfaceSubscriptArray(:,2), aleuroneSurfaceSubscriptArray(:,3)] = ...
%     ind2sub(volumeSize, voxelIndexList);
% 
% aleuroneSurfaceSubscriptArray = aleuroneSurfaceSubscriptArray - grainCenter;
% 
% aleuroneSurfaceSubscriptArray = aleuroneSurfaceSubscriptArray*transform2Vertical;
% 
% aleuroneSurfaceSubscriptArray = aleuroneSurfaceSubscriptArray*transform2Up;



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