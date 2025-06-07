% -------------------------------------------------------------------------
% Script: make_figureS1.m
%
% Description:
% This script generates Figure S1 for the manuscript, comparing the
% readiness potential (RP) in two datasets: original (Park et al., 2020)
% and newly collected data. The comparison is done for the 50% of trials
% closest to the prefered phase of respiration against the remaining 50%.
% The figure contains:
%   - Two readiness potential plots
%
% Data Sources:
%   - Preprocessed EEG and Resp .mat files in data/.../processed folders.
%
% Requirements:
%   - circ_stat toolbox for circular median
%   - Precomputed data in data/.../processed folder (run
%   preprocess_rawdata.m first)
%
% Output:
%   - A composite figure with 1x2 subplots saved as an SVG file suitable 
%     for publication (figure_S1.svg)
% -------------------------------------------------------------------------

clear
close all

%% Load and sort the RPs by participant's prefered phase in both datasets
% % Original data % %

files_original      = dir('data\original\processed\data_EEG_Suj*_E.mat');
RPs_original        = nan(2, 3072, 52);

for participant = 1:length(files_original)
    % Load data
    load([files_original(participant).folder '\' sprintf('data_EEG_Suj%02d_E.mat', participant)])
    load([files_original(participant).folder '\' sprintf('Resp_phase_Suj%02d_E.mat', participant)])
    
    % Eval EEG and Resp
    dat                     = eval(sprintf('data_EEG_Suj%02d_E', participant));
    resp                    = eval(sprintf('Resp_phase_Suj%02d_E', participant));

    % Bin according to Resp phase at time of button press
    Resp_phase_press        = resp(:,2048);
    
    % Get the median value (circular)
    median_value            = circ_median(Resp_phase_press)-pi;

    % Get the length of 25% of the data
    quarter_length          = round(length(Resp_phase_press)/4);

    % Sort the data by respiratory phase
    [sorted_phases, sort_I] = sort(Resp_phase_press);
    
    % Find the median
    pos_median              = find(sorted_phases==median_value);
    % If not exact median find nearest position
    if isempty(pos_median)
        pos_median          = find(sorted_phases>median_value,1,'first');
    end

    % Find indexes for position of 25% before to 25% after the median
    pos_1                   = mod(pos_median - quarter_length - 1, length(Resp_phase_press)) + 1;
    pos_2                   = mod(pos_median + quarter_length - 1, length(Resp_phase_press)) + 1;

    if pos_1<pos_2
        grp_1_I             = sort_I(pos_1:pos_2);
        grp_2_I             = setdiff(sort_I, grp_1_I);
    elseif pos_2<pos_1
        grp_1_I             = [sort_I(pos_1:end)' sort_I(1:pos_2)']';
        grp_2_I             = setdiff(sort_I, grp_1_I);
    else
        error('The start position equals the end position')
    end
    
    % Store result in a data cube
    data_cube               = nan(1,3072,length(sort_I));
    for i = 1:length(sort_I)
        data_cube(:,:,i)    = dat.trial{1,i};
    end

    RPs_original(1,:,participant)   = squeeze(mean(data_cube(1,:,grp_1_I),3)); 
    RPs_original(2,:,participant)   = squeeze(mean(data_cube(1,:,grp_2_I),3));

end

% Save the data for JASP in CSV file
original_rp_diff    = [squeeze(mean(RPs_original(1,(512*2):2048,:),2)) squeeze(mean(RPs_original(2,(512*2):2048,:),2))];
header              = {'phaseA', 'phaseB'};
combinedData        = [header; num2cell(original_rp_diff)];
writecell(combinedData, 'results\jasp_data\original_rp_diff_jasp.csv')

% % New data % %

files_new           = dir('data\new\processed\P*_data_cube.mat');
RPs_new             = nan(2, 2561, 17);

for participant = 1:length(files_new)
    % Load data
    load([files_new(participant).folder '\' sprintf('P%02d_data_cube.mat', participant)])

    % Eval EEG and Resp
    dat                     = squeeze(data_cube(1,:,:));
    resp                    = squeeze(data_cube(2,:,:));

    % Bin according to Resp phase at time of button press
    Resp_phase_press        = resp(2048,:);

    % Get the median value (circular)
    median_value            = circ_median(Resp_phase_press')-pi;

    % Get the length of 25% of the data
    quarter_length          = round(length(Resp_phase_press)/4);

    % Sort the data by respiratory phase
    [sorted_phases, sort_I] = sort(Resp_phase_press);

    % Find the median
    pos_median              = find(sorted_phases==median_value);
    % If not exact median find nearest position
    if isempty(pos_median)
        pos_median          = find(sorted_phases>median_value,1,'first');
    end

    % Find indexes for position of 25% before to 25% after the median
    pos_1                   = mod(pos_median - quarter_length - 1, length(Resp_phase_press)) + 1;
    pos_2                   = mod(pos_median + quarter_length - 1, length(Resp_phase_press)) + 1;

    if pos_1<pos_2
        grp_1_I             = sort_I(pos_1:pos_2);
        grp_2_I             = setdiff(sort_I, grp_1_I);
    elseif pos_2<pos_1
        grp_1_I             = [sort_I(pos_1:end) sort_I(1:pos_2)];
        grp_2_I             = setdiff(sort_I, grp_1_I);
    else
        error('The start position equals the end position')
    end

    RPs_new(1,:,participant)      = squeeze(mean(data_cube(1,:,grp_1_I),3)); 
    RPs_new(2,:,participant)      = squeeze(mean(data_cube(1,:,grp_2_I),3));

end

% Save the data for JASP in CSV file
new_rp_diff         = [squeeze(mean(RPs_new(1,(512*2):2048,:),2)) squeeze(mean(RPs_new(2,(512*2):2048,:),2))];
header              = {'phaseA', 'phaseB'};
combinedData        = [header; num2cell(new_rp_diff)];
writecell(combinedData, 'results\jasp_data\new_rp_diff_jasp.csv')

%% Draw the figure
close all

time                = linspace(-4,1,2560);

fig = figure('Position',[-1548          46        1380         946],'Color','white');

subplot(2,2,[1 2])
upper_grp1_original = squeeze(mean(RPs_original(1,1:2560,:),3)) + squeeze(std(RPs_original(1,1:2560,:),[],3))./sqrt(52);
lower_grp1_original = squeeze(mean(RPs_original(1,1:2560,:),3)) - squeeze(std(RPs_original(1,1:2560,:),[],3))./sqrt(52);

upper_grp2_original = squeeze(mean(RPs_original(2,1:2560,:),3)) + squeeze(std(RPs_original(2,1:2560,:),[],3))./sqrt(52);
lower_grp2_original = squeeze(mean(RPs_original(2,1:2560,:),3)) - squeeze(std(RPs_original(2,1:2560,:),[],3))./sqrt(52);

plot(time, squeeze(mean(RPs_original(1,1:2560,:),3)),'LineWidth',3.5, 'Color',[0 0 0 .2], 'LineStyle','-')
hold on
plot(time, squeeze(mean(RPs_original(2,1:2560,:),3)),'LineWidth',2, 'Color',[0 0 0 1], 'LineStyle','-')

fill([time fliplr(time)], [(lower_grp2_original) fliplr(upper_grp2_original)],[0 0 0],'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.2)
fill([time fliplr(time)], [(lower_grp1_original) fliplr(upper_grp1_original)],[0 0 0],'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.1)

fill([time(1025:2048) fliplr(time(1025:2048))], [(ones(1,1024)*-5) fliplr((ones(1,1024)*2))],[0 0 0],'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.05)


legend('In phase', 'Off phase','Location','southeast','AutoUpdate','off')

ylim([-4.5 2])
yticks([-4.5 -3 -1.5 0 1.5])
xline(0,'--')
ylabel('RP Amplitude (µV)')
text(-.13,0.5,'Park et al. (2020) data','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold','FontSize',30)
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold')
text(1.1,-4.5,'(s)','FontSize',15)
text(1,1,'a','Units','normalized','FontSize',20,'FontWeight','bold')
box off

subplot(2,2,[3 4])
upper_grp1_new = squeeze(mean(RPs_new(1,1:2560,:),3)) + squeeze(std(RPs_new(1,1:2560,:),[],3))./sqrt(17);
lower_grp1_new = squeeze(mean(RPs_new(1,1:2560,:),3)) - squeeze(std(RPs_new(1,1:2560,:),[],3))./sqrt(17);

upper_grp2_new = squeeze(mean(RPs_new(2,1:2560,:),3)) + squeeze(std(RPs_new(2,1:2560,:),[],3))./sqrt(17);
lower_grp2_new = squeeze(mean(RPs_new(2,1:2560,:),3)) - squeeze(std(RPs_new(2,1:2560,:),[],3))./sqrt(17);

plot(time, squeeze(mean(RPs_new(1,1:2560,:),3)),'LineWidth',3.5, 'Color',[0 0 0 .2], 'LineStyle','-')
hold on
plot(time, squeeze(mean(RPs_new(2,1:2560,:),3)),'LineWidth',2, 'Color',[0 0 0 1], 'LineStyle','-')

fill([time fliplr(time)], [(lower_grp1_new) fliplr(upper_grp1_new)],[0 0 0],'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.1)
fill([time fliplr(time)], [(lower_grp2_new) fliplr(upper_grp2_new)],[0 0 0],'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.2)

fill([time(1025:2048) fliplr(time(1025:2048))], [(ones(1,1024)*-5) fliplr((ones(1,1024)*2))],[0 0 0],'LineWidth',2,'EdgeAlpha',0,'FaceAlpha',0.05)

ylim([-4.5 2])
yticks([-4.5 -3 -1.5 0 1.5])
xline(0,'--')
ylabel('RP Amplitude (µV)')
text(-.13,0.5,'New data','Units','normalized','Rotation',90,'HorizontalAlignment','center','FontWeight','bold','FontSize',30)
set(gca,'LineWidth',2,'FontSize',20,'FontWeight','bold')
text(1.1,-4.5,'(s)','FontSize',15)
text(1,1,'b','Units','normalized','FontSize',20,'FontWeight','bold')
box off

print(fig, '-dsvg', '-r300', 'figures\figure_S1.svg', '-painters')