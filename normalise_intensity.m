clear; close all; clc
cd('C:\Users\Admin\Documents\MATLAB\Temp_data')

% Load each data set and save volumes
data =  load('OM16_20_50_2_50_wEndo_distOnFull_voxelFractions_updateSkimming_splitCalcs.mat');

OM1_6_Labels = data.grainVolumeAligned;
OM1_6_Grey = data.greyVolumeAligned;
OM1_6_Index = [data.ALEURONE_INDEX data.ENDOSPERM_INDEX data.GERM_INDEX];

clear data; close all

data =  load('OM17_20_50_2_50_wEndo_distOnFull_voxelFractions_updateSkimming_splitCalcs.mat');

OM1_7_Labels = data.grainVolumeAligned;
OM1_7_Grey = data.greyVolumeAligned;
OM1_7_Index = [data.ALEURONE_INDEX data.ENDOSPERM_INDEX data.GERM_INDEX];

clear data; close all

data =  load('OB6_20_50_2_50_wEndo_distOnFull_voxelFractions_updateSkimming_splitCalcs.mat');

OB6_Labels = data.grainVolumeAligned;
OB6_Grey = data.greyVolumeAligned;
OB6_Index = [data.ALEURONE_INDEX data.ENDOSPERM_INDEX data.GERM_INDEX];

clear data; close all

data =  load('OB7_20_50_2_50_wEndo_distOnFull_voxelFractions_updateSkimming_splitCalcs.mat');

OB7_Labels = data.grainVolumeAligned;
OB7_Grey = data.greyVolumeAligned;
OB7_Index = [data.ALEURONE_INDEX data.ENDOSPERM_INDEX data.GERM_INDEX];

clear data; close all

%% Plot histogram of inensity of all and just germ

figure; subplot(1,3,1); hold on
[n ,x] = hist(OM1_6_Grey(OM1_6_Labels ~= 0), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OM1_7_Grey(OM1_7_Labels ~= 0), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OB6_Grey(OB6_Labels ~= 0), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OB7_Grey(OB7_Labels ~= 0), 1:255);
plot(x, n/sum(n))
title('Full hist')

subplot(1,3,2); hold on
[n ,x] = hist(OM1_6_Grey(OM1_6_Labels == OM1_6_Index(3)), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OM1_7_Grey(OM1_7_Labels == OM1_7_Index(3)), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OB6_Grey(OB6_Labels == OB6_Index(3)), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OB7_Grey(OB7_Labels == OB7_Index(3)), 1:255);
plot(x, n/sum(n))
title('Germ hist')

%% Get means and std
OM1_6_vals = [mean(OM1_6_Grey(OM1_6_Labels ~= 0)) std(single(OM1_6_Grey(OM1_6_Labels ~= 0)))];

OM1_7_vals = [mean(OM1_7_Grey(OM1_7_Labels ~= 0)) std(single(OM1_7_Grey(OM1_7_Labels ~= 0)))];

OB6_vals = [mean(OB6_Grey(OB6_Labels ~= 0)) std(single(OB6_Grey(OB6_Labels ~= 0)))];

OB7_vals = [mean(OB7_Grey(OB7_Labels ~= 0)) std(single(OB7_Grey(OB7_Labels ~= 0)))];

groupMean = mean([OM1_6_vals(1) OM1_7_vals(1) OB6_vals(1) OB7_vals(1)])

groupStd = mean([OM1_6_vals(2) OM1_7_vals(2) OB6_vals(2) OB7_vals(2)])

save('GroupValues.mat', 'OM1_6_vals', 'OM1_7_vals', 'OB6_vals', 'OB7_vals')

%% Test adjustment

OM1_6_Grey2 = uint8((single(OM1_6_Grey)-OM1_6_vals(1))+groupMean);

OM1_7_Grey2 = uint8((single(OM1_7_Grey)-OM1_7_vals(1))+groupMean);

OB6_Grey2 = uint8((single(OB6_Grey)-OB6_vals(1))+groupMean);

OB7_Grey2 = uint8((single(OB7_Grey)-OB7_vals(1))+groupMean);

subplot(1,3,3); hold on
[n ,x] = hist(OM1_6_Grey2(OM1_6_Labels ~= 0), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OM1_7_Grey2(OM1_7_Labels ~= 0), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OB6_Grey2(OB6_Labels ~= 0), 1:255);
plot(x, n/sum(n))

[n ,x] = hist(OB7_Grey2(OB7_Labels ~= 0), 1:255);
plot(x, n/sum(n))
title('Full hist')