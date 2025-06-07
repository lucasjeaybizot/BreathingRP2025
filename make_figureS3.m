% -------------------------------------------------------------------------
% Script: make_figureS3.m
%
% Description:
% This script generates Figure S3 for the manuscript, comparing the
% circular distribution of respiratory phases in newly collected data
% according to 4 separate groupings: W trials, M trials, first 75 trials
% and realigned to EMG. 
% 
% The figure contains:
%   - Four circular histograms
%
% Data Sources:
%   - Raw EEG and Resp .mat files in data/new/raw folder.
%
% Requirements:
%   - four functions GetEMGOnly.m, GetWOnly.m, GetMOnly.m and Get75Only.m
%   in the src\ folder
%
% Output:
%   - A composite figure with 2x2 subplots saved as an SVG file suitable 
%     for publication (figure_S3.svg)
% -------------------------------------------------------------------------

clear
close all

%% Re-preprocess
addpath('src\')

if ~isfile('results\figure_S3\resp_EMG.mat')
    resp_EMG = GetEMG();
else
    load('results\figure_S3\resp_EMG.mat')
end

if ~isfile('results\figure_S3\resp_W.mat')
    resp_W = GetW();
else
    load('results\figure_S3\resp_W.mat')
end

if ~isfile('results\figure_S3\resp_M.mat')
    resp_M = GetM();
else
    load('results\figure_S3\resp_M.mat')
end

if ~isfile('results\figure_S3\resp_75.mat')
    resp_75 = Get75();
else
    load('results\figure_S3\resp_75.mat')
end

%% Draw the figure

fig = figure('Name','Figure 1', 'Position', [65    50   848   946], 'Color','white');

% Plot W trials only
subplot(2,2,1)
total_phases        = [resp_W.resp_phases];
individual_phases   = arrayfun(@(x) circ_mean(x.resp_phases'), resp_W);

h = polarhistogram(total_phases,14,'FaceColor','k','FaceAlpha',0.2); hold on
r = h.Parent.RLim(2);
polarplot(individual_phases, (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k');
polarplot(circ_mean(individual_phases'), (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

text(0.5,1.2,'W trials','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(-0.1,1.1,'a','Units','normalized','FontSize',12,'FontWeight','bold')

% Plot M trials only
subplot(2,2,2)
total_phases        = [resp_M.resp_phases];
individual_phases   = arrayfun(@(x) circ_mean(x.resp_phases'), resp_M);

h = polarhistogram(total_phases,14,'FaceColor','k','FaceAlpha',0.2); hold on
r = h.Parent.RLim(2);
polarplot(individual_phases, (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k');
polarplot(circ_mean(individual_phases'), (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

text(0.5,1.2,'M trials','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(-0.1,1.1,'b','Units','normalized','FontSize',12,'FontWeight','bold')

% Plot first 75 trials only
subplot(2,2,3)
total_phases        = [resp_75.resp_phases];
individual_phases   = arrayfun(@(x) circ_mean(x.resp_phases'), resp_75);

h = polarhistogram(total_phases,14,'FaceColor','k','FaceAlpha',0.2); hold on
r = h.Parent.RLim(2);
polarplot(individual_phases, (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k');
polarplot(circ_mean(individual_phases'), (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

text(0.5,1.2,'First 75 trials','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(-0.1,1.1,'c','Units','normalized','FontSize',12,'FontWeight','bold')

% Plot aligned to EMG onset
subplot(2,2,4)
total_phases        = [resp_EMG.resp_phases];
individual_phases   = arrayfun(@(x) circ_mean(x.resp_phases'), resp_EMG);

h = polarhistogram(total_phases,14,'FaceColor','k','FaceAlpha',0.2); hold on
r = h.Parent.RLim(2);
polarplot(individual_phases, (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k');
polarplot(circ_mean(individual_phases'), (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

text(0.5,1.2,'EMG aligned','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(1.2,0,'N=17','Units','normalized','HorizontalAlignment','center','FontSize',12)
text(-0.1,1.1,'d','Units','normalized','FontSize',12,'FontWeight','bold')

print(fig, '-dsvg', '-r300', 'figures\figure_S3.svg', '-painters')

