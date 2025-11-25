% -------------------------------------------------------------------------
% Script: make_figure1.m
%
% Description:
% This script generates Figure 1 for the manuscript, visualizing the
% relationship between respiration phase and EEG readiness potential (RP)
% amplitude across three datasets: original (Park et al., 2020), newly
% collected, and simulated data. The figure contains:
%   - Rose plots showing the distribution of respiration phases at button
%     press and individual mean directions.
%   - Bar plots of RP amplitude binned by respiration phase, using both
%     original and corrected binning methods. (1b is an exact reproduction
%     of the original figure).
%
% Data Sources:
%   - Results from both original and corrected binning analyses, loaded 
%     from precomputed .mat files in results folder.
%
% Requirements:
%   - circ_stat toolbox for circular statistics
%   - Precomputed data in results folder (run binning_corrected.m and
%   binning_original.m first)
%
% Output:
%   - A composite figure with 3x3 subplots saved as an SVG file suitable 
%     for publication (figure_1.svg)
% -------------------------------------------------------------------------

clear
close all

%% Load processed data for both original and corrected binning methods
load("results\corrected_binning\originaldata_correctedbinning.mat")
load("results\corrected_binning\newdata_correctedbinning.mat")
load("results\corrected_binning\simulateddata_correctedbinning.mat")

load("results\original_binning\originaldata_originalbinning.mat")
load("results\original_binning\newdata_originalbinning.mat")
load("results\original_binning\simulateddata_originalbinning.mat")

%% Initialize figure layout
fig = figure('Name','Figure 1', 'Position', [-1779          50        1380         946], 'Color','white');
% Position for second screen (Used for final figure generation): [-1779          50        1380         946]
% Position for main monitor:      [338          25        1380         946]
% Position for remote second screen: [1681         131        1920
% 955]

%% Rose plots
subplot(3,3,1) % Original data rose plot
total_phases        = vertcat(originaldata_correctedbinning.Resp_phase_press);
individual_phases   = arrayfun(@(x) circ_mean(x.Resp_phase_press), originaldata_correctedbinning);

h = polarhistogram(total_phases,14,'FaceColor','k','FaceAlpha',0.2); hold on
r = h.Parent.RLim(2);
polarplot(individual_phases, (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k');
polarplot(circ_mean(individual_phases'), (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

text(0.5,1.2,'Behavioral coupling','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','normal','FontSize',15)
text(-.3,0.5,'Park et al., 2020','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','normal','FontSize',15)
text(1.2,0,'N=52','Units','normalized','HorizontalAlignment','center','FontSize',12)
text(1.1,1,sprintf('p<%.3f',round(originaldata_originalbinning(1).pval_omni_total,4)),'Units','normalized', 'BackgroundColor', 'w','VerticalAlignment','top','EdgeColor','k')

text(-0.1, 1, 'a', 'Units','normalized', 'FontWeight','normal','FontSize',12)

subplot(3,3,4) % New data rose plot
total_phases        = horzcat(newdata_correctedbinning.Resp_phase_press);
individual_phases   = arrayfun(@(x) circ_mean(x.Resp_phase_press'), newdata_correctedbinning);

h = polarhistogram(total_phases,14,'FaceColor','k','FaceAlpha',0.2); hold on
r = h.Parent.RLim(2);
polarplot(individual_phases, (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k');
polarplot(circ_mean(individual_phases'), (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

text(-.3,0.5,'New data','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','normal','FontSize',15)
text(1.2,0,'N=17','Units','normalized','HorizontalAlignment','center','FontSize',12)
text(1.1,1,sprintf('p=%.3f',round(newdata_originalbinning(1).pval_omni_total,4)),'Units','normalized', 'BackgroundColor', 'w','VerticalAlignment','top','EdgeColor','k')

text(-0.1, 1, 'd', 'Units','normalized', 'FontWeight','normal','FontSize',12)

subplot(3,3,7) % Simulated data rose plot
total_phases        = horzcat(simulateddata_correctedbinning.Resp_phase_press);
individual_phases   = arrayfun(@(x) circ_mean(x.Resp_phase_press'), simulateddata_correctedbinning);

h = polarhistogram(total_phases,6,'FaceColor','k','FaceAlpha',0.2); hold on
r = h.Parent.RLim(2);
polarplot(individual_phases, (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'None', 'MarkerEdgeColor', 'k');
polarplot(circ_mean(individual_phases'), (r)*ones(size(individual_phases)), 'o','MarkerSize',7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

text(-.3,0.5,'Simulated data','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','normal','FontSize',15)
text(1.2,0,'N=50','Units','normalized','HorizontalAlignment','center','FontSize',12)

text(-0.1, 1, 'g', 'Units','normalized', 'FontWeight','normal','FontSize',12)

%% Prep for bar plots
sine_time = linspace(0,7,1000);
sine_wave = sin(linspace((pi/2)-pi/3,(5*pi/2)+pi/3,1000));
sine_wave = sin(linspace((pi/2)-pi/6,(5*pi/2)+pi/6,1000));
%% Original binning bar plots
subplot(3,3,2) % Original data bar plot
mean_amps = reshape(cell2mat([arrayfun(@(x) nanmean(x.mean_amp_sumnorm,1), originaldata_originalbinning, 'UniformOutput', false)]), [], length(originaldata_originalbinning))';
mean_amps = mean_amps(:, [4:6 1:3]);

b = bar(mean(mean_amps,1),'k','FaceAlpha',0.2);
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
scatter(x_bar, mean_amps','k')
% ylim([-0.4 0.8]) % Force ylim to match Park et al. 2020 for easier comparison
y_lim = get(gca,'ylim');
y_lim      = get(gca,'ylim');
range_y    = y_lim(2) - y_lim(1);

bottom     = y_lim(1) + 0.10 * range_y;   % 10% padding bottom
top        = y_lim(1) + 0.90 * range_y;   % 10% padding top

amplitude  = (top - bottom) / 2;
offset     = bottom + amplitude;

sine_to_plot = offset + amplitude * sine_wave;   % centered sine wave

plot(sine_time, sine_to_plot, '--k');

text(0.05,0.965,sprintf('MI: p<%.3f',round(originaldata_originalbinning(1).pval_MI_total,4)),'Units','normalized', 'BackgroundColor', 'w','VerticalAlignment','top','EdgeColor','k')

text(0.5,1.2,'RP coupling: Original binning','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','normal','FontSize',15)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
text(-.15,0.5,'Normalized RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

text(-0.21, 1, 'b', 'Units','normalized', 'FontWeight','normal','FontSize',12)
subplot(3,3,5) % New data bar plot
mean_amps = reshape(cell2mat([arrayfun(@(x) nanmean(x.mean_amp_sumnorm,1), newdata_originalbinning, 'UniformOutput', false)]), [], length(newdata_originalbinning))';
mean_amps = mean_amps(:, [4:6 1:3]);

b = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
scatter(x_bar, mean_amps','k')
y_lim = get(gca,'ylim');
y_lim      = get(gca,'ylim');
range_y    = y_lim(2) - y_lim(1);

bottom     = y_lim(1) + 0.10 * range_y;   % 10% padding bottom
top        = y_lim(1) + 0.90 * range_y;   % 10% padding top

amplitude  = (top - bottom) / 2;
offset     = bottom + amplitude;

sine_to_plot = offset + amplitude * sine_wave;   % centered sine wave

plot(sine_time, sine_to_plot, '--k');

text(0.05,0.965,sprintf('MI: p<%.3f',round(newdata_originalbinning(1).pval_MI_total,4)),'Units','normalized', 'BackgroundColor', 'w','VerticalAlignment','top','EdgeColor','k')

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
text(-.15,0.5,'Normalized RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

text(-0.21, 1, 'e', 'Units','normalized', 'FontWeight','normal','FontSize',12)

subplot(3,3,8) % Simulated data bar plot
mean_amps = reshape(cell2mat([arrayfun(@(x) nanmean(x.mean_amp_sumnorm,1), simulateddata_originalbinning, 'UniformOutput', false)]), [], length(simulateddata_originalbinning))';
mean_amps = mean_amps(:, [4:6 1:3]);

b = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
scatter(x_bar, mean_amps','k')
y_lim = get(gca,'ylim');
y_lim      = get(gca,'ylim');
range_y    = y_lim(2) - y_lim(1);

bottom     = y_lim(1) + 0.10 * range_y;   % 10% padding bottom
top        = y_lim(1) + 0.90 * range_y;   % 10% padding top

amplitude  = (top - bottom) / 2;
offset     = bottom + amplitude;

sine_to_plot = offset + amplitude * sine_wave;   % centered sine wave

plot(sine_time, sine_to_plot, '--k');


xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
text(-.15,0.5,'Normalized RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

text(-0.21, 1, 'h', 'Units','normalized', 'FontWeight','normal','FontSize',12)

%% Corrected binning bar plots
subplot(3,3,3) % Original data bar plot
mean_amps = reshape([originaldata_correctedbinning.mean_amp], [], length(originaldata_correctedbinning))';
mean_amps = mean_amps(:, [4:6 1:3]);

b = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
lowsem =  mean(mean_amps,1)-std(mean_amps,[],1)/sqrt(size(mean_amps,1));
scatter(x_bar, mean_amps','k')
y_lim      = get(gca,'ylim');
range_y    = y_lim(2) - y_lim(1);

bottom     = y_lim(1) + 0.10 * range_y;   % 10% padding bottom
top        = y_lim(1) + 0.90 * range_y;   % 10% padding top

amplitude  = (top - bottom) / 2;
offset     = bottom + amplitude;

sine_to_plot = offset + amplitude * sine_wave;   % centered sine wave

plot(sine_time, sine_to_plot, '--k');

text(0.05,0.965,sprintf('MI: p=%.3f',round(originaldata_correctedbinning(1).pval_MI_total,4)),'Units','normalized', 'BackgroundColor', 'w','VerticalAlignment','top','EdgeColor','k')

text(0.5,1.2,'RP coupling: Corrected binning','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','normal','FontSize',15)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
% ylim([min(sine_to_plot)-(.1*range(sine_to_plot)) max(sine_to_plot)+(.1*range(sine_to_plot))])

annotation('arrow', [0.92 0.92], [50/60 55/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
annotation('arrow', [0.92 0.92], [48/60 43/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
text(1.1, 0.6, 'Inspiration', 'Rotation',90, 'Units', 'normalized')
text(1.1, 0.1, 'Expiration', 'Rotation',90, 'Units', 'normalized')

text(-.15,0.5,'RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

text(-0.1, 1, 'c', 'Units','normalized', 'FontWeight','normal','FontSize',12)

subplot(3,3,6); % New data bar plot
mean_amps = reshape([newdata_correctedbinning.mean_amp], [], length(newdata_correctedbinning))';
mean_amps = mean_amps(:, [4:6 1:3]);

b = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
lowsem =  mean(mean_amps,1)-std(mean_amps,[],1)/sqrt(size(mean_amps,1));
scatter(x_bar, mean_amps', 'k')
y_lim      = get(gca,'ylim');
range_y    = y_lim(2) - y_lim(1);

bottom     = y_lim(1) + 0.10 * range_y;   % 10% padding bottom
top        = y_lim(1) + 0.90 * range_y;   % 10% padding top

amplitude  = (top - bottom) / 2;
offset     = bottom + amplitude;

sine_to_plot = offset + amplitude * sine_wave;   % centered sine wave

plot(sine_time, sine_to_plot, '--k');

text(0.05,0.965,sprintf('MI: p=%.3f',round(newdata_correctedbinning(1).pval_MI_total,4)),'Units','normalized', 'BackgroundColor', 'w','VerticalAlignment','top','EdgeColor','k')

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
ylim([min(sine_to_plot)-(.1*range(sine_to_plot)) max(sine_to_plot)+(.1*range(sine_to_plot))])

annotation('arrow', [0.92 0.92], [50/60 55/60]-18/60, 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
annotation('arrow', [0.92 0.92], [48/60 43/60]-18/60, 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
text(1.1, 0.6, 'Inspiration', 'Rotation',90, 'Units', 'normalized')
text(1.1, 0.1, 'Expiration', 'Rotation',90, 'Units', 'normalized')

text(-.15,0.5,'RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

text(-0.1, 1, 'f', 'Units','normalized', 'FontWeight','normal','FontSize',12)

subplot(3,3,9) % Simulated data bar plot
mean_amps = reshape([simulateddata_correctedbinning.mean_amp], [], length(simulateddata_correctedbinning))';
mean_amps = mean_amps(:, [4:6 1:3]);

b = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
lowsem =  mean(mean_amps,1)-std(mean_amps,[],1)/sqrt(size(mean_amps,1));
scatter(x_bar, mean_amps', 'k')
y_lim      = get(gca,'ylim');
set(gca,'ylim',ylim*1.1)
y_lim      = get(gca,'ylim');
range_y    = y_lim(2) - y_lim(1);

bottom     = y_lim(1) + 0.10 * range_y;   % 10% padding bottom
top        = y_lim(1) + 0.90 * range_y;   % 10% padding top

amplitude  = (top - bottom) / 2;
offset     = bottom + amplitude;

sine_to_plot = offset + amplitude * sine_wave;   % centered sine wave

plot(sine_time, sine_to_plot, '--k');

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
ylim([min(sine_to_plot)-(.1*range(sine_to_plot)) max(sine_to_plot)+(.1*range(sine_to_plot))])

annotation('arrow', [0.92 0.92], [50/60 55/60]-36/60, 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
annotation('arrow', [0.92 0.92], [48/60 43/60]-36/60, 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
text(1.1, 0.6, 'Inspiration', 'Rotation',90, 'Units', 'normalized')
text(1.1, 0.1, 'Expiration', 'Rotation',90, 'Units', 'normalized')

text(-.15,0.5,'RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

text(-0.1, 1, 'i', 'Units','normalized', 'FontWeight','normal','FontSize',12)

print(fig, '-dsvg', '-r300', 'figures\figure_1_scatter.svg', '-painters')