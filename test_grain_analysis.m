% Test load and plot of grain stack.
clc; clear; close all

targetDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Grain LU/Data/Om_1_7_test/Labels';

% Loading tiff stack takes a while.
grainVolume = loadtiffstack(targetDirectory, 1);

volumeSize = size(grainVolume);

aleuroneIndex = 1; % Material index.



%% Take mesh around entire grain then limit to aleurone exterior.

% Reduce patch size to save mesh construction time.
% If this value is below ~0.05, holes may appear in mesh if grain touches volume border.
resizeRatio = 0.2;

smallGrainVolume = imresize3(grainVolume, resizeRatio);

smallGrainVolume = imbinarize(smallGrainVolume, 0);

% Do a strong morphological open to seperate along crease.
%smallGrainVolume = imopen(smallGrainVolume, strel('disk', 3));

%%% Should test for presence of multiple objects here.

% Create isosurface on downsampled volume.
fullGrainSurface = isosurface(smallGrainVolume, 0.5);

%Simplify patch to speed processing time later. 
patchReductionRatio = 0.01;

fullGrainSurface = reducepatch(fullGrainSurface,patchReductionRatio);

% Englarge surface back to original size.
fullGrainSurface.vertices = round(fullGrainSurface.vertices/resizeRatio);

% Note that surface vertices are flipped relative to image subscripts.
fullGrainSurface.vertices(:,[1 2]) = fullGrainSurface.vertices(:,[2 1]);

% Truncate to lower bound of volume.
fullGrainSurface.vertices(fullGrainSurface.vertices < 1) = 1;
        
% Truncate to upper bound.
fullGrainSurface.vertices(fullGrainSurface.vertices(:,1) > volumeSize(1),1) = volumeSize(1);

fullGrainSurface.vertices(fullGrainSurface.vertices(:,2) > volumeSize(2),2) = volumeSize(2);

fullGrainSurface.vertices(fullGrainSurface.vertices(:,3) > volumeSize(3),3) = volumeSize(3);

vertexNo = size(fullGrainSurface.vertices,1);

faceNo = size(fullGrainSurface.faces,1);

% Plost test figures;
figure; axis equal; hold on; axis off; set(gca, 'Clipping', 'off')

patch(fullGrainSurface);

plot3(fullGrainSurface.vertices(:,1), fullGrainSurface.vertices(:,2), fullGrainSurface.vertices(:,3), 'ro')

title('Check there are no holes in the surface')



%% Get voxels on exterior of orignal grain by overlap to exterior.
%%% Make function for this and replace in loop.
grainExterior = ~grainVolume;

grainExterior = imdilate(grainExterior, strel('disk', 1));

grainExterior = grainExterior & grainVolume;

% Get surface indices, then convert to subscripts.
grainSurfaceVoxelIndexList = find(grainExterior);

indexNo = length(grainSurfaceVoxelIndexList);

grainSurfaceSubscriptArray = zeros(indexNo, 3);

[grainSurfaceSubscriptArray(:,1), grainSurfaceSubscriptArray(:,2), grainSurfaceSubscriptArray(:,3)] = ...
    ind2sub(volumeSize, grainSurfaceVoxelIndexList);



%% Snap each surface vertices onto the surface voxels and test if aleurone.
vertexToRemove = zeros(vertexNo,1);

facesToRemove = zeros(faceNo,1);

%%%voxelRange = 1/resizeRatio;

% If voxel range is to small some subvolumes may be empty.
%%%if voxelRange <= 20; voxelRange = voxelRange * 1.5; end

voxelRange = 25;

for iVertex = 1:vertexNo
    surfaceVertex = fullGrainSurface.vertices(iVertex,:);
    
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
        
        %%% Replace in volume test with function and add above
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

            fullGrainSurface.vertices(iVertex,:) = surfaceVertex;
            
        else
            
           error('No surface voxels in subvolume') 
           
        end
    end
    
    
    
    % Test if vertex is part of aleurone and mark to remove if not.
    voxelMaterialIndex = grainVolume(surfaceVertex(1), surfaceVertex(2), surfaceVertex(3));
    
    if aleuroneIndex ~= voxelMaterialIndex
        vertexToRemove(iVertex) = 1;
        
        % Also mark faces to remove.
        vertexInFace = find(fullGrainSurface.faces == iVertex);
        
        [tempFaceIndex, ~] = ind2sub([faceNo, 3], vertexInFace);
        
        facesToRemove(tempFaceIndex) = 1;
    end
end

%%% Make function for vertex/face removal and use below
% Remove faces, vertices, and update numbers.
facesToRemove = find(facesToRemove);

fullGrainSurface.faces(facesToRemove, :) = [];

vertexToRemove = find(vertexToRemove);

% Need to keep original vertex positions to relink faces
initialVertexPositions = 1:vertexNo;

fullGrainSurface.vertices(vertexToRemove, :) = [];

initialVertexPositions(vertexToRemove) = [];

vertexNo = size(fullGrainSurface.vertices,1);

faceNo = size(fullGrainSurface.faces,1);

% Relink remaining faces.
for iVertex = 1:vertexNo
   
    fullGrainSurface.faces(fullGrainSurface.faces == initialVertexPositions(iVertex)) = iVertex;
    
end

% Plot test figures.
figure; axis equal; hold on; axis off; set(gca, 'Clipping', 'off')

patch(fullGrainSurface); 

plot3(fullGrainSurface.vertices(:,1), fullGrainSurface.vertices(:,2), fullGrainSurface.vertices(:,3), 'ro')



%% Test geodesic embedding for unwrapping.
% Problem can occur if there are unlinked vertices, check by running fast marching once.
tempD = perform_fast_marching_mesh(fullGrainSurface.vertices, fullGrainSurface.faces, 50);

%%% If small number of unreachable faces, should be able to find all by testing
%%% one point randomly. Need to confirm.

% Remove points with inf distance.
if sum(isinf(tempD))
    vertexToRemove = find(isinf(tempD));

    facesToRemove = zeros(faceNo,1);

    for iVertex = 1:length(vertexToRemove)
        vertexInFace = find(fullGrainSurface.faces == vertexToRemove(iVertex));

        [tempFaceIndex, ~] = ind2sub([faceNo, 3], vertexInFace);

        facesToRemove(tempFaceIndex) = 1;
    end

    facesToRemove = find(facesToRemove);

    fullGrainSurface.faces(facesToRemove, :) = [];

    % Need to keep original vertex positions to relink faces
    initialVertexPositions = 1:vertexNo;

    fullGrainSurface.vertices(vertexToRemove, :) = [];

    initialVertexPositions(vertexToRemove) = [];

    vertexNo = size(fullGrainSurface.vertices,1);

    faceNo = size(fullGrainSurface.faces,1);

    % Relink remaining faces.
    for iVertex = 1:vertexNo

        fullGrainSurface.faces(fullGrainSurface.faces == initialVertexPositions(iVertex)) = iVertex;

    end
end

% Get distance map for all.
D = zeros(vertexNo);

for iVertex = 1:vertexNo
    
    D(:,iVertex) = perform_fast_marching_mesh(fullGrainSurface.vertices, fullGrainSurface.faces, iVertex);
    
end

sum(isinf(D(:)))

D = (D + D')/2;

J = eye(vertexNo) - ones(vertexNo)/vertexNo;
W = -J*(D.^2)*J;

%%% Infs are sometimes produced in loop, need to resolve.

[U,S] = eig(W);
S = diag(S);
[S,I] = sort(S,'descend'); U = U(:,I);

% To plot.
vertexF = U(:,1:2)' .* repmat(sqrt(S(1:2)), [1 vertexNo]);

icenter = 50;
irotate = 50;

vertexF = vertexF - repmat(vertexF(:,icenter), [1 vertexNo]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];

figure;
plot_mesh(vertexF,fullGrainSurface.faces);
       
% To do: 
%        - Then apply deformation field to volumes
%        - Look down into volume to get thickness of aleurone
%        - Could also integrate intensity to see if this varies