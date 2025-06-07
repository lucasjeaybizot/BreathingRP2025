% -------------------------------------------------------------------------
% Script: make_figureS2.m
%
% Description:
% This script generates Figure S2 for the manuscript, comparing the
% circular distribution of respiratory phases in two datasets: original
% (Park et al., 2020) and newly collected data. The comparison is done for
% the distribution of respiratory phases at time of button press, versus
% across randomly shuffled distributions of respiratory phases. The figure
% contains:
%   - Two circular histograms (one per dataset)
%
% Data Sources:
%   - Preprocessed EEG and Resp .mat files in data/.../processed folders.
%
% Requirements:
%   - Precomputed data in data/.../processed folder (run
%   preprocess_rawdata.m first)
%
% Output:
%   - A composite figure with 1x2 subplots saved as an SVG file suitable 
%     for publication (figure_S2.svg)
% -------------------------------------------------------------------------

clear
close all

%% Get distributions from both New and Original datasets
% % New data % %
files_new                   = dir("data\new\processed\P*_resp_segment.mat");
resp_phase_new              = [];
resp_perm_new               = [];

for participant = 1:length(files_new)
    % Load the respiration data
    load(fullfile(files_new(participant).folder, files_new(participant).name))

    % Rotate respiration to account for flip when saving .eeg
    resp_segment_tp         = resp_segment;
    resp_segment(resp_segment_tp>0) = resp_segment_tp(resp_segment_tp>0) - pi;
    resp_segment(resp_segment_tp<0) = resp_segment_tp(resp_segment_tp<0) + pi;

    % Append respiration phases at time of press
    resp_phase_new          = [resp_phase_new resp_segment(press_positions)];

    % Get permutted respiration phases
    tp_perm                 = nan(1000, length(press_positions));
    for i = 1:1000
        tp_seg                  = circshift(resp_segment, randi(length(resp_segment)));
        tp_perm(i,:)            = tp_seg(press_positions);
    end
    resp_perm_new           = [resp_perm_new tp_perm(:)'];
end

% % Original data % %
files_original              = dir("data\original\processed\Resp_phase_Suj*_E.mat");
resp_phase_original         = [];
resp_perm_original          = [];

for participant = 1:length(files_original)
    % Load the respiration data
    load(fullfile(files_original(participant).folder, files_original(participant).name))

    % Due to data available being epoched, concatenate to allow shuffling
    resp                = eval(sprintf('Resp_phase_Suj%02d_E', participant));
    cat_resp            = reshape(resp.', 1, []);

    % Extract and append respiratory phase at time of press
    resp_phase_original = [resp_phase_original cat_resp(2048:3073:end)];

    % Get permutted respiration phases
    tp_perm                 = nan(1000, length(cat_resp(2048:3073:end)));
    for i = 1:1000
        tp_seg                  = circshift(cat_resp, randi(length(cat_resp)));
        tp_perm(i,:)            = tp_seg(2048:3073:end);
    end
    resp_perm_original           = [resp_perm_original tp_perm(:)'];
end

%% Plot the figure
close all

fig = figure('Position', [680   307   910   571], 'Color', 'white');

% First subplot
pax1 = polaraxes; % Create a unique polar axes for the first subplot
hold(pax1, 'on');

% Plot the first polar histogram with transparency
h1 = polarhistogram(pax1, resp_phase_new, 'Normalization', 'probability', ...
                    'NumBins', 14, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the second polar histogram with different color and transparency
h2 = polarhistogram(pax1, resp_perm_new, 'Normalization', 'probability', ...
                    'NumBins', 14, 'FaceColor', 'black', 'FaceAlpha', 0.6, 'EdgeColor', 'none');

% Draw a reference circle at a specific radius
radius = 0.0637; % Set the radius for the reference circle
theta = linspace(0, 2*pi, 100); % Full circle
polarplot(pax1, theta, radius * ones(size(theta)), 'k--', 'LineWidth', 1.5);

% Adjust radial limits based on the histogram values
rmax1 = max([h1.Values, h2.Values]);
rlim(pax1, [0 .1]);

pax1.Position = [0.55 0.1 0.35 0.8]; % Adjust as needed

text(-1.37,1.15,'a','Units','normalized','FontSize',12,'FontWeight','bold')

hold(pax1, 'off');

% Second subplot
pax2 = polaraxes; % Create a unique polar axes for the second subplot
hold(pax2, 'on');

% Plot the first polar histogram with transparency
h3 = polarhistogram(pax2, resp_phase_original, 'Normalization', 'probability', ...
                    'NumBins', 14, 'FaceColor', 'black', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the second polar histogram with different color and transparency
h4 = polarhistogram(pax2, resp_perm_original, 'Normalization', 'probability', ...
                    'NumBins', 14, 'FaceColor', 'black', 'FaceAlpha', 0.6, 'EdgeColor', 'none');

% Draw a reference circle at a specific radius
radius2 = min([h4.Values]); % Reference circle radius based on h4 minimum
polarplot(pax2, theta, radius2 * ones(size(theta)), 'k--', 'LineWidth', 1.5);

% Adjust radial limits based on the histogram values
rmax2 = max([h3.Values, h4.Values]);
rlim(pax2, [0 .1]);

pax2.Position = [0.1 0.1 0.35 0.8];
leg1 = legend(pax2,'actual distribution','shuffled distribution');
set(leg1, 'Position', [0.6, 0.05, 0.3, 0.05]);

title(pax2, 'Orignal data')
title(pax1, 'New data')
text(1.2,1.15,'b','Units','normalized','FontSize',12,'FontWeight','bold')
hold(pax2, 'off');

print(fig, '-dsvg', '-r300', 'figures\figure_S2.svg', '-painters')