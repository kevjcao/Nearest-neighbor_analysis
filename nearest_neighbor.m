%% Calculates the average for the nearest neighbor of every cell and displays histogram of average distances between cells 
% Input: FIJI/ImageJ files from Measurements that report the ROI number, 
% ROI area, X coordinate, and Y coordinate in columns 1, 2, 3, and 4, 
% respectively. 

% Output: Histogram of intercell distances for ROIs, heatmap of Euclidean
% distances, and the average distance between the nearest neighboring cell.

% NOTE: this code requires the Image Processing Toolbox

%% Clear workspace
clear;
close all;

%% Import file
[FileName, FilePath] = uigetfile( ...
       {'*.xlsx','New Excel (*.xlsx)'; ...
       '*.xls','Old Excel (*.xls)';
       '*.txt', 'Text file (*.txt)'}, ...
        'Pick a file', ...
        'MultiSelect', 'on');
File = fullfile(FilePath, FileName);
[num, txt, raw] = xlsread(File);

fprintf('%s \n\n', FileName);                                               % Display file name / experiment details
tic
% clearvars txt FilePath FileName raw

%% Parameters
BW_dist = 10;                                                               % in microns, bin size for all euclidean distances b/n ROIs
BW_nn = 5;                                                                  % in microns, bin size for histogram of the nearest neighbors
scale = 1.1;                                                                % pixels / micron (from micrometer)
overlap = 0;                                                                % 1 if overlap of cells is allowed

%% Convert raw ROIs from pixels to microns
% Datatype: column 1 = ROI number, 2 = area of ROI, 3 = X coord, 4 = Y coord
% Create matrix that contains only X,Y coordinates for pdist function
raw_coord(:,1) = num(:,3);
raw_coord(:,2) = num(:,4);

raw_coord_m = raw_coord / scale;                                            % convert raw_coord from pixels to microns

%% Convert measured area to radius
soma_array = num(:,2) / scale^2;
rad_array = sqrt(soma_array(:,1) / pi);
Rad.mean = mean(rad_array.');
Rad.std = std(rad_array.');
fprintf('Average radius of light responsive cells = %.2f +/- %.2f microns \n\n', Rad.mean, Rad.std);

%% Plot ROIs
scatter_labels = num2str(num(:,1));

coord_plot = figure(1);
hold on
% plot((num(:,3) / scale), (num(:,4) / scale), 'b*');
% text((num(:,3) / scale), (num(:,4) / scale), scatter_labels);
viscircles(raw_coord_m, rad_array);
text(raw_coord_m(:,1), raw_coord_m(:,2), scatter_labels);
grid on;
axis equal
xlim([0 (256 / scale)]);
ylim([0 (256 / scale)]);
set(gca, 'YDir', 'reverse');
xlabel('X coordinate, microns');
ylabel('Y coordinate, microns');
set(gca, 'TickDir', 'out');
box on;
hold off

%% Calculate the Euclidean distances between all ROIS
% get Euclidean distances b/n every point
eucDis = pdist(raw_coord_m, 'euclidean');

eucDis_histo = figure(2);
hold on
histogram(eucDis, 'BinWidth', BW_dist);
xlabel('Euclidean distance between cells, microns');
ylabel('Counts');
set(gca, 'TickDir', 'out');
hold off

% plot heatmap of the measured Euclidean distances
eucDis_squareform = squareform(eucDis);
eucDisSQform = figure(3);
heatmap(eucDis_squareform);


%% Find nearest neighboring light-responsive cell
% minimum Euclidean distances
eucDis_min = eucDis_squareform;
eucDis_min(eucDis_min == 0) = NaN;
nearestNeighb = min(eucDis_min);

% mean of all nearest neighbors
NNeighb.mean = mean(nearestNeighb);
NNeighb.std = std(nearestNeighb);
fprintf('Average distance between nearest light-responsive cells = %.2f +/- %.2f microns \n\n', NNeighb.mean, NNeighb.std);

% nearest neighbor histo
[N_nn, Edges_nn] = histcounts(nearestNeighb, 'BinWidth', BW_nn);

nn_histo = figure(4);
hold on
histogram(nearestNeighb, 'BinWidth', BW_nn);
xl = xline(NNeighb.mean, '--', {'Average'});
xl.LineWidth = 2;
xlabel('Distance of nearest light-responsive cell, microns');
ylabel('Counts');
set(gca, 'TickDir', 'out');
hold off

%% Tiling check
% randomize x & y coordinates for each ROI

shuffle_px(:,1) = randi([0 256], 1, length(num(:,2)));
shuffle_px(:,2) = randi([0 256], 1, length(num(:,3)));
shuffle_coord = shuffle_px / scale;                                             

shuffle_plot = figure(5);                                                   % first pass shuffle
hold on
viscircles(shuffle_coord, rad_array);
text(shuffle_coord(:,1), shuffle_coord(:,2), scatter_labels);
axis equal
xlim([0 (256 / scale)]);
ylim([0 (256 / scale)]);
set(gca, 'YDir', 'reverse');
xlabel('X coordinate, microns');
ylabel('Y coordinate, microns');
set(gca, 'TickDir', 'out');
hold off


% measure euclidean distances between nearest neighbor of shuffled ROIs
shuf_eucDis = pdist(shuffle_coord, 'euclidean');
shuf_min = squareform(shuf_eucDis);
shuf_min(shuf_min == 0) = NaN;
[M, I] = min(shuf_min);
check_shuf_dist(:,1) = I.';                                                 % cell index for nearest neighbor
check_shuf_dist(:,2) = M.';                                                 % euclidean distance between nearest neighbor
check_shuf_dist(:,3) = rad_array + rad_array(check_shuf_dist(:,1));         % minimum distance between two cells without overlap

% reshuffle if overlap is present (e.g. euclidean distance between nearest 
% neighbor < sum of raddi for nearest neighbors)

if overlap == 0                                                             % = 1 if analysis conditions allow cell overlap
    for i = 1:length(check_shuf_dist)
        while check_shuf_dist(i,2) < check_shuf_dist(i,3)
            shuffle_coord(i,1) = randi([0 256], 1) / scale;
            shuffle_coord(i,2) = randi([0 256], 1) / scale;
            shuf_eucDis = pdist(shuffle_coord, 'euclidean');
            shuf_min = squareform(shuf_eucDis);
            shuf_min(shuf_min == 0) = NaN;
            [M, I] = min(shuf_min);
            check_shuf_dist(:,1) = I.';                                         % cell index for nearest neighbor
            check_shuf_dist(:,2) = M.';                                         % euclidean distance between nearest neighbor
            check_shuf_dist(:,3) = rad_array + rad_array(check_shuf_dist(:,1)); % minimum distance between two cells without overlap
            if check_shuf_dist(i,2) > check_shuf_dist(i,3)
                i = i + 1;
            end
        end
    end
else
end

figure(6)
heatmap(shuf_min);

shuf_nn = min(shuf_min);

% mean of all nearest neighbors
NNshuf.mean = mean(shuf_nn);
NNshuf.std = std(shuf_nn);
fprintf('Average distance between nearest shuffled light-responsive cells = %.2f +/- %.2f microns \n\n', NNshuf.mean, NNshuf.std);

% nearest neighbor histo
[N_nn_shuf, Edges_nn_shuf] = histcounts(shuf_nn, 'BinWidth', BW_nn);

nn_histo = figure(7);
hold on
histogram(shuf_nn, 'BinWidth', BW_nn);
xl = xline(NNshuf.mean, '--', {'Average'});
xl.LineWidth = 2;
xlabel('Distance of nearest shuffled light-responsive cell, microns');
ylabel('Counts');
set(gca, 'TickDir', 'out');
hold off


% fprintf('No overlapping cells in shuffle \n\n');
shuffle_plot = figure(8);
hold on
viscircles(shuffle_coord, rad_array);
text(shuffle_coord(:,1), shuffle_coord(:,2), scatter_labels);
axis equal
grid on;
xlim([0 (256 / scale)]);
ylim([0 (256 / scale)]);
set(gca, 'YDir', 'reverse');
xlabel('X coordinate, microns');
ylabel('Y coordinate, microns');
set(gca, 'TickDir', 'out');
box on;
hold off

toc

