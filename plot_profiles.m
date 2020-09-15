%Load and plot combined profiles on one map

clear; close all; clc
cd('C:\Users\Admin\Documents\MATLAB\Temp_data\PlotProfiles')

OB6_data = load('OB6_profiles.mat');

OB7_data = load('OB7_profiles.mat');

OM16_data = load('OM16_profiles.mat');

OM17_data = load('OM17_profiles.mat');

%%

%Plot vertical
figure

subplot(1,2,1); hold on
plot(OB6_data.verticalThicknessMean, OB6_data.xLabelVert, 'b');
plot(OB7_data.verticalThicknessMean, OB7_data.xLabelVert, 'r');
plot(OM16_data.verticalThicknessMean, OM16_data.xLabelVert, '--','color', [0 0 0.75]);
plot(OM17_data.verticalThicknessMean, OM17_data.xLabelVert, '--','color', [0.75 0 0]);

xlabel('Average aleurone thickness (um)')
ylabel('Distance from base (mm)')
xlim([0 60]); ylim([0 12])

subplot(1,2,2); hold on
plot(OB6_data.verticalIntensityMean, OB6_data.xLabelVert);
plot(OB7_data.verticalIntensityMean, OB7_data.xLabelVert, 'r');
plot(OM16_data.verticalIntensityMean, OM16_data.xLabelVert, '--','color', [0 0 0.75]);
plot(OM17_data.verticalIntensityMean, OM17_data.xLabelVert, '--','color', [0.75 0 0]);

xlabel('Average aleurone intensity (AU)')
title('Vertical profiles')
ylabel('Distance from base (mm)')
xlim([0 200]); ylim([0 12])

% subplot(1,3,3); hold on
% plot(OB6_data.verticalCount*OB6_data.areaPerPoint, OB6_data.xLabelVert);
% plot(OB7_data.verticalCount*OB7_data.areaPerPoint, OB7_data.xLabelVert, 'r');
% plot(OM16_data.verticalCount*OM16_data.areaPerPoint, OM16_data.xLabelVert, '--','color', [0 0 0.75]);
% plot(OM17_data.verticalCount*OM17_data.areaPerPoint, OM17_data.xLabelVert, '--','color', [0.75 0 0]);
% 
% xlabel('Average aleurone area (mm2)')
% xlim([0 10]); ylim([0 12])

%Plot hoirzontal
figure

subplot(1,2,1); hold on
plot(OB6_data.xLabelHoriz, OB6_data.horizontalThicknessMean, 'b');
plot(OB7_data.xLabelHoriz, OB7_data.horizontalThicknessMean, 'r');
plot(OM16_data.xLabelHoriz, OM16_data.horizontalThicknessMean, '--','color', [0 0 0.75]);
plot(OM17_data.xLabelHoriz, OM17_data.horizontalThicknessMean, '--','color', [0.75 0 0]);

ylabel('Average aleurone thickness (um)')
ylim([0 60]); xlim([-5 5])
xlabel('Distance from centreline (mm)')

subplot(1,2,2); hold on
plot(OB6_data.xLabelHoriz, OB6_data.horizontalIntensityMean);
plot(OB7_data.xLabelHoriz, OB7_data.horizontalIntensityMean, 'r');
plot(OM16_data.xLabelHoriz, OM16_data.horizontalIntensityMean, '--','color', [0 0 0.75]);
plot(OM17_data.xLabelHoriz, OM17_data.horizontalIntensityMean, '--','color', [0.75 0 0]);

ylabel('Average aleurone intensity (AU)')
title('Horizontal profiles')
xlabel('Distance from centreline (mm)')
ylim([0 200]); xlim([-5 5])

legend('OB6', 'OB7', 'OM1-6', 'OM1-7')

% subplot(1,3,3); hold on
% plot(OB6_data.xLabelHoriz, OB6_data.horizontalCount*OB6_data.areaPerPoint);
% plot(OB7_data.xLabelHoriz, OB7_data.horizontalCount*OB7_data.areaPerPoint, 'r');
% plot(OM16_data.xLabelHoriz, OM16_data.horizontalCount*OM16_data.areaPerPoint, '--','color', [0 0 0.75]);
% plot(OM17_data.xLabelHoriz, OM17_data.horizontalCount*OM17_data.areaPerPoint, '--','color', [0.75 0 0]);
% 
% ylabel('Average aleurone area (mm2)')
% xlim([-5 5]); ylim([0 5])
