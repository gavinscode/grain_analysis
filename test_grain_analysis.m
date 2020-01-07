% Test load and plot of grain stack.
clc; clear; close all

targetDirectory = '/Users/gavintaylor/Documents/Company/Client Projects/Grain LU/Data/Om_1_7_test/Labels';

% Loading tiff stack takes a while.
grainVolume = loadtiffstack(targetDirectory, 1);

%% Take mesh around entire grain then limit to aleurone exterior.

% Resize smaller to save processing time.
smallGrainVolume = imresize3(grainVolume, 0.1);

% Create isosurface on downsampled volume.
fullGrainSurface = isosurface(smallGrainVolume, 0.5);

% To do: - Go back and get indices of voxels on surface of grain
%        - Snap vertices of surface onto those voxels
%        - Remove vertices that are not aleurone from faces list

%        - Then look at unwrapping remaining (not closed mesh)
%        - Try techniques from numerical tours
%        - Then apply deformation field to volumeä
%        - Look down into volume to get thickness of aleurone
%        - Could also integrate intensity to see if this varies
%% Plot test figures.

figure;

subplot(2,1,1); patch(fullGrainSurface)

subplot(2,1,2); imshow(smallGrainVolume(:,:,50)*100)