% -------------------------------------------------------------------------
% Script: make_figureS5.m
%
% Description:
% This script generates Figure S5 for the manuscript, comparing original
% and corrected binning bar plots when using sum division normalization
% versus not.
% 
% The figure contains:
%   - Four bar plots
%
% Data Sources:
%   - original binned original data in results/original_binning folder.
%   - corrected binned original data in results/corrected_binning folder.
%
% Requirements:
%   - N/A
%
% Output:
%   - A composite figure with 2x2 subplots saved as an SVG file suitable 
%     for publication (figure_S5.svg)
% -------------------------------------------------------------------------

clear
close all

%% Load the data
load('results\original_binning\originaldata_originalbinning.mat')
load('results\corrected_binning\originaldata_correctedbinning.mat')

%% Draw the figure

fig                 = figure('Name','Figure 1', 'Position', [-1910          50        1005         844], 'Color','white');

sine_time           = linspace(0,7,1000);
sine_wave           = sin(linspace((pi/2)-pi/3,(5*pi/2)+pi/3,1000));

% % original binning - sum normalization (same as Fig 3B in Park et al. 2020) % %
subplot(2,2,1)

mean_amps           = reshape(cell2mat([arrayfun(@(x) nanmean(x.mean_amp_sumnorm,1), originaldata_originalbinning, 'UniformOutput', false)]), [], length(originaldata_originalbinning))';
mean_amps           = mean_amps(:, [4:6 1:3]);

b                   = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar               = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
ylim([-0.4 0.8]) % Force ylim to match Park et al. 2020 for easier comparison
y_lim = get(gca,'ylim');
sine_to_plot        = (range(y_lim)*.2)/2+y_lim(1)+((sine_wave+1)/2)*range(y_lim)*0.8;
plot(sine_time,  sine_to_plot, '--k');

text(0.5,1.1,'Sum normalization','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
text(-.15,0.5,'Normalized RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)
text(-.3,0.5,'Original binning','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(-0.2, 1, 'a', 'FontSize',12,'FontWeight','bold','Units','normalized')

% % corrected binning - sum normalization % % 
subplot(2,2,3)

mean_amps           = reshape([originaldata_correctedbinning.mean_amp_sumnorm], [], length(originaldata_correctedbinning))';
mean_amps           = mean_amps(:, [4:6 1:3]);

b                   = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar               = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
y_lim               = get(gca,'ylim');
sine_to_plot        = (range(y_lim)*.2)/2+y_lim(1)+((sine_wave+1)/2)*range(y_lim)*0.8;
plot(sine_time,  sine_to_plot, '--k');
text(-.15,0.5,'Normalized RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);

text(-.3,0.5,'Corrected binning','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(-0.2, 1, 'c', 'FontSize',12,'FontWeight','bold','Units','normalized')

% % original binning - no normalization % %
subplot(2,2,2)

mean_amps           = reshape([originaldata_originalbinning.mean_amp], [], length(originaldata_originalbinning))';
mean_amps           = mean_amps(:, [4:6 1:3]);

b                   = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar               = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
lowsem              =  mean(mean_amps,1)-std(mean_amps,[],1)/sqrt(size(mean_amps,1));
sine_to_plot        = -1.1*max(abs(mean(lowsem,1)))+((sine_wave+1)/2) * 1.2 * abs(max(abs(mean(lowsem,1))));
plot(sine_time,  sine_to_plot, '--k');

text(0.5,1.1,'No normalization','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(-.15,0.5,'RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
ylim([min(sine_to_plot)-(.1*range(sine_to_plot)) max(sine_to_plot)+(.1*range(sine_to_plot))])

text(-0.1, 1, 'b', 'FontSize',12,'FontWeight','bold','Units','normalized')
annotation('arrow', [0.92 0.92], [47/60 54/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
annotation('arrow', [0.92 0.92], [43/60 36/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
text(1.1, 0.65, 'Inspiration', 'Rotation',90, 'Units', 'normalized')
text(1.1, 0.15, 'Expiration', 'Rotation',90, 'Units', 'normalized')

% % corrected binning - no normalization (same as figure 1c) % %
subplot(2,2,4)
mean_amps       = reshape([originaldata_correctedbinning.mean_amp], [], length(originaldata_correctedbinning))';
mean_amps       = mean_amps(:, [4:6 1:3]);

b               = bar(mean(mean_amps,1),'k','FaceAlpha',0.2); 
hold on
x_bar           = b.XData ;
errorbar(x_bar, mean(mean_amps,1),std(mean_amps,[],1)/sqrt(size(mean_amps,1)), 'k', 'linestyle', 'none','LineWidth',1);
lowsem          =  mean(mean_amps,1)-std(mean_amps,[],1)/sqrt(size(mean_amps,1));
sine_to_plot    = -1.1*max(abs(mean(lowsem,1)))+((sine_wave+1)/2) * 1.2 * abs(max(abs(mean(lowsem,1))));
plot(sine_time,  sine_to_plot, '--k');

text(-.15,0.5,'RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
ylim([min(sine_to_plot)-(.1*range(sine_to_plot)) max(sine_to_plot)+(.1*range(sine_to_plot))])
text(-0.1, 1, 'd', 'FontSize',12,'FontWeight','bold','Units','normalized')
annotation('arrow', [0.92 0.92], [19/60 26/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
annotation('arrow', [0.92 0.92], [15/60 8/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
text(1.1, 0.65, 'Inspiration', 'Rotation',90, 'Units', 'normalized')
text(1.1, 0.15, 'Expiration', 'Rotation',90, 'Units', 'normalized')

print(fig, '-dsvg', '-r300', 'figures\figure_S5.svg', '-painters')