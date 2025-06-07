% -------------------------------------------------------------------------
% Script: make_figureS4.m
%
% Description:
% This script generates Figure S4 for the manuscript, illustrating the
% impact of removing a single outlier (single trial, or single participant). 
% 
% The figure contains:
%   - Three bar plots (Figure 3B from Park et al. 2020)
%
% Data Sources:
%   - original binned original data in results/original_binning folder.
%
% Requirements:
%   - N/A
%
% Output:
%   - A composite figure with 1x3 subplots saved as an SVG file suitable 
%     for publication (figure_S4.svg)
% -------------------------------------------------------------------------

clear
close all

%% Remove outliers (single most extreme values)
% Load data
load('results\original_binning\originaldata_originalbinning.mat')

[~,part_I_sn]           = max(arrayfun(@(x) max(abs(x.mean_amp_sumnorm), [],'all'), originaldata_originalbinning));  % gets the participant with the max value
[~, trl_I_sn]           = max(abs(originaldata_originalbinning(part_I_sn).mean_amp_sumnorm),[],1);                   % gets the max in each bin
trl_rm_I_sn             = mode(trl_I_sn);                                                                            % gets the trial with the most maximums

% Recompute mean amplitude without most extreme value (single trial)
mean_amps_sn            = reshape(cell2mat([arrayfun(@(x) nanmean(x.mean_amp_sumnorm,1), originaldata_originalbinning, 'UniformOutput', false)]), [], length(originaldata_originalbinning))';
mean_amps_sn            = mean_amps_sn(:, [4:6 1:3]);

originaldata_originalbinning(part_I_sn).mean_amp_sumnorm(trl_rm_I_sn,:) = nan(6,1);

mean_amps_outlier_sn    = reshape(cell2mat([arrayfun(@(x) nanmean(x.mean_amp_sumnorm,1), originaldata_originalbinning, 'UniformOutput', false)]), [], length(originaldata_originalbinning))';
mean_amps_outlier_sn    = mean_amps_outlier_sn(:, [4:6 1:3]);

% Find & remove single most extreme participant (single participant)
[~,I_part]              = max(range(mean_amps_sn,2));
mean_amps_outlier_sn_part = mean_amps_sn(setdiff(1:52,I_part),:);

%% Draw the figure
fig = figure('Name','Figure 1', 'Position', [-1779          50        (1380*(3/2))*(8.5/10)         (946/2)*(8.5/10)], 'Color','white');

% % Subplot - left most % %
subplot(1,3,1)
sine_time = linspace(0,7,1000);
sine_wave = sin(linspace((pi/2)-pi/6,(5*pi/2)+pi/6,1000));

b = bar(mean(mean_amps_sn,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps_sn,1),std(mean_amps_sn,[],1)/sqrt(size(mean_amps_sn,1)), 'k', 'linestyle', 'none','LineWidth',1);
ylim([-0.4 0.8]) % Force ylim to match Park et al. 2020 for easier comparison
y_lim = get(gca,'ylim');
sine_to_plot    = (range(y_lim)*.2)/2+y_lim(1)+((sine_wave+1)/2)*range(y_lim)*0.8;
plot(sine_time,  sine_to_plot, '--k');

text(0.5,1.05,'No outlier removal','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
text(-.15,0.5,'Normalized RP amplitude (μV)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontSize',12)
text(-.25,0.5,'Figure 3 B (Park et al. 2020)','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)
text(-0.1, 1, 'a', 'FontSize',12,'FontWeight','bold','Units','normalized')

% % Subplot - center % %
subplot(1,3,2)
sine_time = linspace(0,7,1000);
sine_wave = sin(linspace((pi/2)-pi/6,(5*pi/2)+pi/6,1000));

b = bar(mean(mean_amps_outlier_sn,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps_outlier_sn,1),std(mean_amps_outlier_sn,[],1)/sqrt(size(mean_amps_outlier_sn,1)), 'k', 'linestyle', 'none','LineWidth',1);
ylim([-0.4 0.8]) % Force ylim to match Park et al. 2020 for easier comparison
y_lim = get(gca,'ylim');
sine_to_plot    = (range(y_lim)*.2)/2+y_lim(1)+((sine_wave+1)/2)*range(y_lim)*0.8;
plot(sine_time,  sine_to_plot, '--k');

text(0.5,1.05,'Single trial outlier removed','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
text(-0.1, 1, 'b', 'FontSize',12,'FontWeight','bold','Units','normalized')

% % Subplot - right most % %
subplot(1,3,3)
sine_time = linspace(0,7,1000);
sine_wave = sin(linspace((pi/2)-pi/6,(5*pi/2)+pi/6,1000));

b = bar(mean(mean_amps_outlier_sn_part,1),'k','FaceAlpha',0.2); 
hold on
x_bar = b.XData ;
errorbar(x_bar, mean(mean_amps_outlier_sn_part,1),std(mean_amps_outlier_sn_part,[],1)/sqrt(size(mean_amps_outlier_sn_part,1)), 'k', 'linestyle', 'none','LineWidth',1);
ylim([-0.4 0.8]) % Force ylim to match Park et al. 2020 for easier comparison
y_lim = get(gca,'ylim');
sine_to_plot    = (range(y_lim)*.2)/2+y_lim(1)+((sine_wave+1)/2)*range(y_lim)*0.8;
plot(sine_time,  sine_to_plot, '--k');

text(0.5,1.05,'Single participant outlier removal','Units','normalized','Rotation',0,'HorizontalAlignment','center','FontWeight','bold','FontSize',15)

xticks([1:6]);
xticklabels({'0~60','','','','', '300~360'});
xtickangle(0);
text(-0.1, 1, 'c', 'FontSize',12,'FontWeight','bold','Units','normalized')
annotation('arrow', [0.92 0.92], [32/60 50/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
annotation('arrow', [0.92 0.92], [28/60 10/60], 'Units', 'normalized','LineStyle', '--','HeadStyle','vback3');
text(1.1, 0.6, 'Inspiration', 'Rotation',90, 'Units', 'normalized')
text(1.1, 0.2, 'Expiration', 'Rotation',90, 'Units', 'normalized')

print(fig, '-dsvg', '-r300', 'figures\figure_S4.svg', '-painters')
